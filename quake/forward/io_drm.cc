/* Independent single-file DRM sampling. The three-stage replay workflow in
 * drm.cc and its mesh/refinement rules are deliberately separate. */
#include <fstream>
#include <cstring>
#include <memory>
#include <sys/stat.h>
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
struct Recorder {
    int print_rate = 0, local_count = 0;
    int32_t* device_nodes = NULL;
    int32_t* device_modes = NULL;
    double* device_coords = NULL;
    double* device_values = NULL;
    drm_output::FlatWriter writer;
};
std::vector<std::unique_ptr<Recorder> > recorders;

bool same_file(const std::string& a, const std::string& b) {
    if (a.empty() || b.empty()) return false;
    if (a == b) return true;
    struct stat first, second;
    return stat(a.c_str(), &first) == 0 && stat(b.c_str(), &second) == 0 &&
           first.st_dev == second.st_dev && first.st_ino == second.st_ino;
}

void cuda_check(cudaError_t code, const char* operation) {
    if (code != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(code));
}

void upload(void** device, const void* host, size_t bytes) {
    cuda_check(cudaMalloc(device, bytes), "DRM metadata allocation");
    cuda_check(cudaMemcpy(*device, host, bytes, cudaMemcpyHostToDevice), "DRM metadata copy");
}

// Setup one output at a time: only the largest global point array is retained
// on the host. GPU metadata and MPI views persist independently for each file.
void setup_recording(int32_t rank, std::vector<drm_output::Point>& points,
                     const std::string& path, int print_rate, int batch_size,
                     double surface_shift, int topography,
                     double domain_x, double domain_y, double domain_z,
                     mysolver_t* solver, double delta_t) {
    char output_path[256] = "";
    int config[3] = {0, print_rate, batch_size};
    if (rank == 0) {
        if (path.empty() || path.size() > 240)
            throw std::runtime_error("DRM output paths must contain 1 to 240 bytes");
        std::strcpy(output_path, path.c_str());
        for (size_t i = 0; i < points.size(); ++i) {
            drm_output::Point& p = points[i];
            if (!std::isfinite(p[0]) || !std::isfinite(p[1]) ||
                p[0] < 0 || p[0] >= domain_x || p[1] < 0 || p[1] >= domain_y)
                throw std::runtime_error(path + ": point " + std::to_string(i) +
                    " lies outside the horizontal domain; box offsets locate its top center");
            p[2] += surface_shift;
            // Match station/plane surface-relative depth conventions.
            if (topography) p[2] += point_elevation(p[0], p[1]);
            if (!std::isfinite(p[2]) || p[2] < 0 || p[2] >= domain_z)
                throw std::runtime_error(path + ": point " + std::to_string(i) + " lies outside the vertical domain");
        }
        config[0] = static_cast<int>(points.size());
        printf("DRM output: %s; %d samples, XYZ displacement in metres; dt_out=%.17g s, batch=%d snapshots\n",
               output_path, config[0], print_rate * delta_t, batch_size);
    }
    drm_output::mpi_check(MPI_Bcast(config, 3, MPI_INT, 0, comm_solver), "DRM config broadcast");
    recorders.emplace_back(new Recorder);
    Recorder& recorder = *recorders.back();
    recorder.print_rate = config[1];
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
    recorder.local_count = static_cast<int>(ids.size());
    if (recorder.local_count) {
        if (!solver || !solver->gpuData || !solver->gpuData->tm1Device)
            throw std::runtime_error("DRM flat output requires device displacement state");
        upload(reinterpret_cast<void**>(&recorder.device_nodes), nodes.data(), nodes.size() * sizeof(int32_t));
        upload(reinterpret_cast<void**>(&recorder.device_modes), modes.data(), modes.size() * sizeof(int32_t));
        upload(reinterpret_cast<void**>(&recorder.device_coords), coords.data(), coords.size() * sizeof(double));
        cuda_check(cudaMalloc(reinterpret_cast<void**>(&recorder.device_values),
                              static_cast<size_t>(recorder.local_count) * 3 * sizeof(double)), "DRM result allocation");
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
    recorder.writer.open(comm_solver, output_path, config[0], ids, config[2]);
}

void setup(int32_t rank, const char* parameters, double surface_shift,
           int topography, double* corner_lon, double* corner_lat,
           double domain_x, double domain_y, double domain_z,
           UTMZone_t* zone, mysolver_t* solver, double delta_t) {
    int counts[2] = {0, 0}; // generated boxes, optional CSV recorder
    int box_rate = 0, csv_rate = 0, batch_size = 16;
    char directory[256] = "", csv_output[256] = "", csv_input[256] = "";
    std::vector<drm_output::BoxSpec> boxes;
    if (rank == 0) {
        FILE* fp = hu_fopen(parameters, "r");
        if (hu_config_get_int_opt(fp, "number_output_drm_boxes", &counts[0]) < 0 || counts[0] < 0)
            throw std::runtime_error("number_output_drm_boxes must be a nonnegative integer");
        char* path = csv_output;
        size_t length = sizeof(csv_output);
        counts[1] = hu_config_get_string_opt(fp, "drm_output_file", &path, &length);
        if (counts[1] < 0 || (counts[1] && !csv_output[0]))
            throw std::runtime_error("Invalid drm_output_file");
        if (counts[0] || counts[1]) {
            if (hu_config_get_int_opt(fp, "drm_output_batch_size", &batch_size) < 0 || batch_size <= 0)
                throw std::runtime_error("drm_output_batch_size must be a positive integer");
        }
        if (counts[0]) {
            if (hu_config_get_int_req(fp, "output_drm_boxes_print_rate", &box_rate) != 1 || box_rate <= 0)
                throw std::runtime_error("output_drm_boxes_print_rate must be a positive integer");
            path = directory;
            length = sizeof(directory);
            if (hu_config_get_string_req(fp, "output_drm_boxes_directory", &path, &length) != 1 || !directory[0])
                throw std::runtime_error("output_drm_boxes_directory is required");
            std::ifstream input(parameters);
            boxes = drm_output::read_boxes(input, "output_drm_boxes", counts[0]);
        }
        if (counts[1]) {
            if (hu_config_get_int_req(fp, "drm_print_rate", &csv_rate) != 1 || csv_rate <= 0)
                throw std::runtime_error("drm_print_rate must be a positive integer");
            path = csv_input;
            length = sizeof(csv_input);
            if (hu_config_get_string_opt(fp, "drm_points_file", &path, &length) != 1 || !csv_input[0])
                throw std::runtime_error("drm_output_file requires drm_points_file; use output_drm_boxes for generated boxes");
        }
        fclose(fp);
    }
    drm_output::mpi_check(MPI_Bcast(counts, 2, MPI_INT, 0, comm_solver), "DRM recorder counts");
    std::vector<std::string> opened_paths;
    const size_t total = static_cast<size_t>(counts[0]) + counts[1];
    for (size_t index = 0; index < total; ++index) {
        std::vector<drm_output::Point> points;
        std::string path;
        const bool is_box = index < static_cast<size_t>(counts[0]);
        if (rank == 0) {
            if (is_box) {
                path = std::string(directory) + "/drmdisplacements." + std::to_string(index);
                drm_output::BoxGeometry geometry = drm_output::box_faces(boxes[index]);
                printf("DRM box %zu: top center=(%.17g,%.17g), exterior=(%.17g,%.17g,%.17g), "
                       "spacing=%.17g m, strike=%.17g deg clockwise from north\n", index,
                       boxes[index][0], boxes[index][1], boxes[index][2], boxes[index][3],
                       boxes[index][4], boxes[index][5], boxes[index][6]);
                for (const drm_output::Face& face : geometry.faces)
                    printf("  %s: rows=%d columns=%d point_offset=%zu\n",
                           face.name, face.rows, face.columns, face.offset);
                points.swap(geometry.points);
            } else {
                path = csv_output;
                std::ifstream input(csv_input);
                if (!input) throw std::runtime_error(std::string("Cannot open DRM points: ") + csv_input);
                points = drm_output::read_points(input);
                for (drm_output::Point& p : points) {
                    vector3D_t xy = compute_domain_coords_linearinterp(
                        p[1], p[0], corner_lon, corner_lat, domain_y, domain_x, zone);
                    p[0] = xy.x[0];
                    p[1] = xy.x[1];
                }
                printf("DRM CSV recorder: preserving input row order\n");
            }
            if (same_file(path, parameters) || same_file(path, csv_input))
                throw std::runtime_error("DRM output would overwrite an input file: " + path);
            for (const std::string& opened : opened_paths)
                if (same_file(path, opened))
                    throw std::runtime_error("DRM recorders must use distinct files: " + path);
            opened_paths.push_back(path);
        }
        setup_recording(rank, points, path, is_box ? box_rate : csv_rate, batch_size,
                        surface_shift, topography, domain_x, domain_y, domain_z, solver, delta_t);
    }
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
    try {
        for (const auto& entry : recorders) {
            Recorder& r = *entry;
            if (step % r.print_rate) continue;
            if (r.local_count) {
                const cudaStream_t stream = solver->streams[CUDA_STREAM_MAIN];
                psolve_launch_plane_interpolate((r.local_count + 255) / 256, 256, 0, stream,
                    r.local_count, r.device_nodes, r.device_coords, r.device_modes,
                    solver->gpuData->tm1Device, r.device_values);
                cuda_check(cudaGetLastError(), "DRM interpolation kernel");
                const size_t bytes = static_cast<size_t>(r.local_count) * 3 * sizeof(double);
                cuda_check(cudaMemcpyAsync(r.writer.next_frame(), r.device_values, bytes,
                                          cudaMemcpyDeviceToHost, stream), "DRM packed result copy");
                cuda_check(cudaStreamSynchronize(stream), "DRM result synchronization");
                solver->gpu_spec->numbytespci += bytes;
            }
            r.writer.finish_frame();
        }
    } catch (const std::exception& error) {
        solver_abort("drm_output_print", NULL, "%s", error.what());
    }
}

void drm_output_close() {
    try {
        for (const auto& entry : recorders) {
            Recorder& r = *entry;
            r.writer.close();
            cuda_check(cudaFree(r.device_nodes), "DRM nodes free");
            cuda_check(cudaFree(r.device_modes), "DRM modes free");
            cuda_check(cudaFree(r.device_coords), "DRM coordinates free");
            cuda_check(cudaFree(r.device_values), "DRM values free");
        }
        recorders.clear();
    } catch (const std::exception& error) {
        solver_abort("drm_output_close", NULL, "%s", error.what());
    }
}
