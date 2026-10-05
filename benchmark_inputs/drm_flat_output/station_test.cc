// Actual production station interpolation/dispatch on an analytic host field.
// Compile with STATION_TEST_MPI to also use the production collective writer.
#include <cassert>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>
#ifdef STATION_TEST_MPI
#include "../../quake/forward/flat_output_mpi.h"
static flat_output::FlatWriter stationBinaryWriter;
#else
struct MemoryWriter {
    std::vector<double> frame, output;
    double* next_frame() { return frame.data(); }
    void finish_frame() { output.insert(output.end(), frame.begin(), frame.end()); }
    void close() {}
} stationBinaryWriter;
#endif
#define __FUNCTION_NAME __func__
const int YES = 1, NO = 0, FEM = 0;
struct vector3D_t { double x[3]; };
struct fvector_t { double f[3]; };
struct tensor_t { double xx, yy, zz, xy, yz, xz; };
struct mysolver_t { fvector_t tm1[8], tm2[8]; } solver;
struct Station {
    int id;
    vector3D_t localcoords;
    int nodestointerpolate[8];
    int pressure_is_defined;
    double acoustic_disp[3], acoustic_vel[3], acoustic_acc[3], prev_disp[3];
    double h, mu, lambda;
    FILE* fpoutputfile;
};
struct TopoStation {
    int TopoStation, nodes_to_interpolate[4];
    double local_coord[3], prev_disp[3];
};
static std::vector<TopoStation> myTopoStations;
static std::vector<std::vector<char> > stationAsciiBuffers;
static int theTopoMethod = 1;
static struct {
    int stationBinary, myNumberOfStations, theNumberOfStations;
    int includeTopography, printStationVelocities, printStationAccelerations;
    int printStationStrains, printStationStresses, theStationsPrintRate;
    double theDeltaT, theDeltaTSquared;
    Station* myStations;
    const char* theStationsDirOut;
} Param;
static struct { mysolver_t* mySolver; } Global = {&solver};
static std::atomic<int> ancillaryCalls(0);
static void Timer_Start(const char*) {}
static void Timer_Stop(const char*) {}
static void solver_station_tensor_kick(int) { ++ancillaryCalls; }
static void solver_sync_pore_pressure_station_values(mysolver_t*) { ++ancillaryCalls; }
static void station_write_pressure_columns(int) { ++ancillaryCalls; }
static const double* solver_station_strain_get(int) { assert(false); return NULL; }
static const double* solver_station_ve_stress_get(int) { assert(false); return NULL; }
static int compute_tetra_strain(int, fvector_t*, double, tensor_t*) { assert(false); return 0; }
static tensor_t point_strain(fvector_t*, double, double, double, double) { assert(false); return {}; }
static tensor_t point_stress(tensor_t, double, double) { assert(false); return {}; }
static void hu_fclose(FILE* fp) { assert(fclose(fp) == 0); }
[[noreturn]] static void solver_abort(const char*, const char*, const char*, ...) {
    throw std::runtime_error("unexpected solver abort");
}
#ifdef _OPENMP
static int station_output_use_omp() { return 1; }
#endif
#include "station_functions.inc"

static double field(int step, int axis, double x, double y, double z) {
    return 100*step + 10*axis + (axis+1)*x + (axis+2)*y + (axis+3)*z;
}
static double expected(int step, int id, int axis) {
    if (id == 2) return 100*step + 3*axis + 0.125; // Reconstructed acoustic state.
    if (id == 1) return field(step, axis, -0.4, -0.6, -0.8); // Cut-cell tet.
    return field(step, axis, 0.25, -0.5, 0.75); // Hex8; id 3 repeats id 0.
}
int main(int argc, char** argv) {
    int rank = 0, ranks = 1;
#ifdef STATION_TEST_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    // Reject independent singleton worlds before any shared output is opened.
    const int expected_ranks = argc == 3 && std::strcmp(argv[2], "1") == 0 ? 1 :
                               argc == 3 && std::strcmp(argv[2], "4") == 0 ? 4 : 0;
    if (!expected_ranks || ranks != expected_ranks) {
        if (rank == 0) {
            if (!expected_ranks)
                fprintf(stderr, "Usage: %s OUTPUT_PATH EXPECTED_RANKS(1|4)\n", argv[0]);
            else
                fprintf(stderr, "MPI world size mismatch: expected %d, got %d; check the MPI compiler and launcher\n",
                        expected_ranks, ranks);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (rank == 0) printf("Station output: verified MPI world size %d\n", ranks);
#else
    (void)argc; (void)argv;
#endif
    const int active = ranks == 1 ? 1 : ranks - 1; // Include an empty rank.
    std::vector<int> ids;
    for (int id = 0; id < 4; ++id) if (id % active == rank) ids.push_back(id);
    std::vector<Station> stations(ids.size());
    myTopoStations.resize(ids.size());
    Param.theNumberOfStations = 4;
    Param.myNumberOfStations = ids.size();
    Param.myStations = stations.empty() ? NULL : stations.data();
    Param.includeTopography = YES;
    Param.theStationsPrintRate = 2;
    Param.theDeltaT = 0.01;
    Param.theDeltaTSquared = 0.0001;
    for (size_t i = 0; i < ids.size(); ++i) {
        Station& st = stations[i];
        st.localcoords = {{0.25, -0.5, 0.75}};
        for (int n = 0; n < 8; ++n) st.nodestointerpolate[n] = n;
        st.pressure_is_defined = ids[i] == 2;
        if (ids[i] == 1) myTopoStations[i] = {1, {0,1,2,4}, {0.2,0.3,0.1}, {0,0,0}};
    }
    // Repeat modes: default ASCII first, binary second. File handles remain
    // NULL in binary mode, exposing any accidental ASCII/ancillary writes.
    for (int binary = NO; binary <= YES; ++binary) {
        assert(Param.stationBinary == NO); // Zero initialization preserves ASCII default.
        Param.stationBinary = binary;
        ancillaryCalls = 0;
        if (binary) {
#ifdef STATION_TEST_MPI
            stationBinaryWriter.open(MPI_COMM_WORLD, argv[1], 4, ids, 1);
#else
            stationBinaryWriter.frame.resize(ids.size()*3);
#endif
        } else {
            for (auto& st : stations) { st.fpoutputfile = tmpfile(); assert(st.fpoutputfile); }
        }
        for (int step = 0; step < 10; ++step) {
            for (int n = 0; n < 8; ++n)
                for (int axis = 0; axis < 3; ++axis)
                    solver.tm1[n].f[axis] = field(step, axis, (n&1)?1:-1, (n&2)?1:-1, (n&4)?1:-1);
            for (size_t i = 0; i < ids.size(); ++i)
                for (int axis = 0; axis < 3; ++axis)
                    stations[i].acoustic_disp[axis] = expected(step, ids[i], axis);
            solver_output_stations(step);
        }
        if (!binary) {
            for (size_t i = 0; i < ids.size(); ++i) {
                rewind(stations[i].fpoutputfile);
                for (int step = 0; step < 10; step += 2) {
                    double time, u[3];
                    assert(fscanf(stations[i].fpoutputfile, "%lf %lf %lf %lf", &time, &u[0], &u[1], &u[2]) == 4);
                    assert(std::abs(time - step*Param.theDeltaT) < 1e-12);
                    for (int axis = 0; axis < 3; ++axis)
                        assert(std::abs(u[axis] - expected(step, ids[i], axis)) < 1e-10);
                }
                double extra;
                assert(fscanf(stations[i].fpoutputfile, "%lf", &extra) == EOF);
            }
        } else assert(ancillaryCalls == 0);
        solver_close_station_outputs();
    }
#ifdef STATION_TEST_MPI
    if (rank == 0) {
        FILE* fp = fopen(argv[1], "rb"); assert(fp);
        for (int step = 0; step < 10; step += 2)
            for (int id = 0; id < 4; ++id)
                for (int axis = 0; axis < 3; ++axis) {
                    double value;
                    assert(fread(&value, sizeof(value), 1, fp) == 1);
                    assert(std::abs(value - expected(step, id, axis)) < 1e-10);
                }
        assert(fgetc(fp) == EOF);
        fclose(fp);
    }
    MPI_Finalize();
#else
    assert(stationBinaryWriter.output.size() == 5*4*3);
    size_t offset = 0;
    for (int step = 0; step < 10; step += 2)
        for (int id = 0; id < 4; ++id)
            for (int axis = 0; axis < 3; ++axis)
                assert(std::abs(stationBinaryWriter.output[offset++] - expected(step, id, axis)) < 1e-10);
#endif
    if (rank == 0) puts("Station ASCII/binary displacement and sampling regression passed");
}
