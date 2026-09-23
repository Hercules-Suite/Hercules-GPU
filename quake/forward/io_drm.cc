/* Independent single-file DRM sampling. The three-stage replay workflow in
 * drm.cc and its mesh/refinement rules are deliberately separate. */
#include <fstream>
#include <cstring>
#include "drm_output_points.h"
#include "drm_output_mpi.h"
#include "io_drm.h"
#include "psolve_gpu_launch.h"
#include "topography.h"
#include "util.h"
#include "quake_util.h"

extern int compute_csi_eta_dzeta(octant_t*, vector3D_t, vector3D_t*, int32_t*);
extern void createDir(const char*, const char*);

namespace {
bool enabled = false;
int print_rate = 0;
int local_count = 0;
int32_t* device_nodes = NULL;
int32_t* device_modes = NULL;
double* device_coords = NULL;
double* device_values = NULL;
drm_output::FlatWriter writer;

void cuda_check(cudaError_t code, const char* operation) {
    if (code != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(code));
}

void upload(void** device, const void* host, size_t bytes) {
    cuda_check(cudaMalloc(device, bytes), "DRM metadata allocation");
    cuda_check(cudaMemcpy(*device, host, bytes, cudaMemcpyHostToDevice), "DRM metadata copy");
}

void setup(int32_t rank, const char* parameters, double surface_shift,
           int topography, double* corner_lon, double* corner_lat,
           double domain_x, double domain_y, double domain_z,
           UTMZone_t* zone, mysolver_t* solver, double delta_t) {
    char output_path[256] = "";
    int config[3] = {0, 0, 16}; // point count, print rate, batch size
    std::vector<drm_output::Point> points;
    if (rank == 0) {
        FILE* fp = hu_fopen(parameters, "r");
        char* path = output_path;
        size_t length = sizeof(output_path);
        const int found = hu_config_get_string_opt(fp, "drm_output_file", &path, &length);
        if (found < 0 || (found && (!output_path[0] || strlen(output_path) > 240)))
            throw std::runtime_error("drm_output_file must be a nonempty path of at most 240 bytes");
        if (found) {
            if (hu_config_get_int_req(fp, "drm_print_rate", &config[1]) != 1 || config[1] <= 0)
                throw std::runtime_error("drm_print_rate must be a positive integer");
            if (hu_config_get_int_opt(fp, "drm_output_batch_size", &config[2]) < 0 || config[2] <= 0)
                throw std::runtime_error("drm_output_batch_size must be a positive integer");
            char points_path[256] = "";
            path = points_path;
            length = sizeof(points_path);
            const int from_file = hu_config_get_string_opt(fp, "drm_points_file", &path, &length);
            if (from_file < 0 || (from_file && !points_path[0]))
                throw std::runtime_error("Invalid drm_points_file path (maximum 255 bytes)");
            if (from_file) {
                std::ifstream input(points_path);
                if (!input) throw std::runtime_error(std::string("Cannot open DRM points: ") + points_path);
                points = drm_output::read_points(input);
                for (size_t i = 0; i < points.size(); ++i) {
                    vector3D_t xy = compute_domain_coords_linearinterp(
                        points[i][1], points[i][0], corner_lon, corner_lat,
                        domain_y, domain_x, zone);
                    points[i][0] = xy.x[0];
                    points[i][1] = xy.x[1];
                }
            } else {
                double boundary[5], spacing, offset_x, offset_y;
                if (parsetext(fp, "drm_edgesize", 'd', &spacing) != 0 ||
                    parsetext(fp, "drm_offset_x", 'd', &offset_x) != 0 ||
                    parsetext(fp, "drm_offset_y", 'd', &offset_y) != 0 ||
                    parsedarray(fp, "drm_boundary", 5, boundary) != 0)
                    throw std::runtime_error("Generated DRM output requires drm_edgesize, drm_offset_x/y and drm_boundary");
                points = drm_output::box_points(boundary, spacing, offset_x, offset_y);
            }
            for (size_t i = 0; i < points.size(); ++i) {
                drm_output::Point& p = points[i];
                if (!std::isfinite(p[0]) || !std::isfinite(p[1]) ||
                    p[0] < 0 || p[0] >= domain_x || p[1] < 0 || p[1] >= domain_y)
                    throw std::runtime_error("DRM point " + std::to_string(i) + " lies outside the horizontal domain");
                p[2] += surface_shift;
                // Match station/plane depth conventions, including terrain.
                if (topography) p[2] += point_elevation(p[0], p[1]);
                if (!std::isfinite(p[2]) || p[2] < 0 || p[2] >= domain_z)
                    throw std::runtime_error("DRM point " + std::to_string(i) + " lies outside the vertical domain");
            }
            config[0] = static_cast<int>(points.size());
            printf("DRM flat output: %s; %d points (%s), XYZ displacement in metres; "
                   "dt_out=%.17g s, batch=%d snapshots\n", output_path, config[0],
                   from_file ? "input row order" : "inner then outer surfaces, X/Y/Z order",
                   config[1] * delta_t, config[2]);
        }
        fclose(fp);
    }
    drm_output::mpi_check(MPI_Bcast(config, 3, MPI_INT, 0, comm_solver), "DRM config broadcast");
    if (!config[0]) return;
    enabled = true;
    print_rate = config[1];
    points.resize(config[0]);
    static_assert(sizeof(drm_output::Point) == 3 * sizeof(double), "Coordinate triples must be packed");
    drm_output::mpi_check(MPI_Bcast(points.data(), config[0] * 3, MPI_DOUBLE, 0, comm_solver),
                          "DRM coordinate broadcast");
    drm_output::mpi_check(MPI_Bcast(output_path, sizeof(output_path), MPI_CHAR, 0, comm_solver),
                          "DRM path broadcast");

    std::vector<int> owners(config[0], 0), ids;
    std::vector<int32_t> nodes, modes;
    std::vector<double> coords;
    for (int i = 0; i < config[0]; ++i) {
        vector3D_t point = {{points[i][0], points[i][1], points[i][2]}};
        octant_t* octant = NULL;
        if (!search_point(point, &octant)) continue;
        edata_t material;
        if (solver_get_octant_material(octant, &material) < 0 ||
            !(material.Vp > 0) || !(material.rho > 0) || !(material.Vs > ACOUSTIC_VS_EPS))
            throw std::runtime_error("DRM point " + std::to_string(i) +
                " is not in a solid element; scalar-acoustic displacement reconstruction is unsupported");
        int32_t stencil[8];
        vector3D_t local;
        if (!compute_csi_eta_dzeta(octant, point, &local, stencil))
            throw std::runtime_error("Cannot interpolate DRM point " + std::to_string(i));
        const bool tetra = stencil[4] == 0 && stencil[5] == 0 && stencil[6] == 0 && stencil[7] == 0;
        double sum = 0;
        for (int axis = 0; axis < 3; ++axis) {
            const double c = local.x[axis];
            if (!std::isfinite(c) || c < (tetra ? -1e-7 : -1.0000001) || c > 1.0000001)
                throw std::runtime_error("Invalid interpolation coordinates at DRM point " + std::to_string(i));
            sum += c;
        }
        if (tetra && sum > 1.0000001)
            throw std::runtime_error("DRM point outside solid topography tetrahedron: " + std::to_string(i));
        owners[i] = 1;
        ids.push_back(i);
        nodes.insert(nodes.end(), stencil, stencil + 8);
        coords.insert(coords.end(), local.x, local.x + 3);
        modes.push_back(tetra ? PLANE_GPU_POINT_TOPO : PLANE_GPU_POINT_TRILINEAR);
    }
    drm_output::mpi_check(MPI_Allreduce(MPI_IN_PLACE, owners.data(), config[0], MPI_INT,
                                       MPI_SUM, comm_solver), "DRM ownership check");
    for (int i = 0; i < config[0]; ++i)
        if (owners[i] != 1)
            throw std::runtime_error("DRM point " + std::to_string(i) + " has " +
                                     std::to_string(owners[i]) + " owners; expected exactly one");
    local_count = static_cast<int>(ids.size());
    if (local_count) {
        if (!solver || !solver->gpuData || !solver->gpuData->tm1Device)
            throw std::runtime_error("DRM flat output requires device displacement state");
        upload(reinterpret_cast<void**>(&device_nodes), nodes.data(), nodes.size() * sizeof(int32_t));
        upload(reinterpret_cast<void**>(&device_modes), modes.data(), modes.size() * sizeof(int32_t));
        upload(reinterpret_cast<void**>(&device_coords), coords.data(), coords.size() * sizeof(double));
        cuda_check(cudaMalloc(reinterpret_cast<void**>(&device_values),
                              static_cast<size_t>(local_count) * 3 * sizeof(double)), "DRM result allocation");
    }
    if (rank == 0) {
        std::string parent(output_path);
        const size_t slash = parent.find_last_of('/');
        if (slash != std::string::npos)
            createDir(parent.substr(0, slash ? slash : 1).c_str(), "drm_output_setup");
        printf("DRM flat backend: GPU plane interpolation + collective MPI-IO on solver ranks; "
               "solid displacement (poroelastic skeleton where applicable); no mesh refinement\n");
    }
    drm_output::mpi_check(MPI_Barrier(comm_solver), "DRM directory creation");
    writer.open(comm_solver, output_path, config[0], ids, config[2]);
}
} // namespace

void drm_output_setup(int32_t rank, const char* parameters, double surface_shift,
                      int topography, double* corner_lon, double* corner_lat,
                      double domain_x, double domain_y, double domain_z,
                      UTMZone_t* zone, mysolver_t* solver, double delta_t) {
    try {
        setup(rank, parameters, surface_shift, topography, corner_lon, corner_lat,
              domain_x, domain_y, domain_z, zone, solver, delta_t);
    } catch (const std::exception& error) {
        solver_abort("drm_output_setup", NULL, "%s", error.what());
    }
}

void drm_output_print(mysolver_t* solver, int32_t step) {
    if (!enabled || step % print_rate) return;
    try {
        if (local_count) {
            const cudaStream_t stream = solver->streams[CUDA_STREAM_MAIN];
            psolve_launch_plane_interpolate((local_count + 255) / 256, 256, 0, stream,
                local_count, device_nodes, device_coords, device_modes,
                solver->gpuData->tm1Device, device_values);
            cuda_check(cudaGetLastError(), "DRM interpolation kernel");
            const size_t bytes = static_cast<size_t>(local_count) * 3 * sizeof(double);
            cuda_check(cudaMemcpyAsync(writer.next_frame(), device_values, bytes,
                                      cudaMemcpyDeviceToHost, stream), "DRM packed result copy");
            cuda_check(cudaStreamSynchronize(stream), "DRM result synchronization");
            solver->gpu_spec->numbytespci += bytes;
        }
        writer.finish_frame();
    } catch (const std::exception& error) {
        solver_abort("drm_output_print", NULL, "%s", error.what());
    }
}

void drm_output_close() {
    if (!enabled) return;
    try {
        writer.close();
        cuda_check(cudaFree(device_nodes), "DRM nodes free");
        cuda_check(cudaFree(device_modes), "DRM modes free");
        cuda_check(cudaFree(device_coords), "DRM coordinates free");
        cuda_check(cudaFree(device_values), "DRM values free");
        enabled = false;
    } catch (const std::exception& error) {
        solver_abort("drm_output_close", NULL, "%s", error.what());
    }
}
