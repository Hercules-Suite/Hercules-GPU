// Real collective I/O regression; no CUDA or solver build required.
#include "../../quake/forward/flat_output_mpi.h"
#include "../../quake/forward/drm_output_points.h"
#include "../../quake/forward/drm_output_coordinates.h"
#include "../../quake/forward/drm_output_materials_mpi.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    // Require an explicit count: an incompatible MPI launcher can silently
    // start multiple singleton worlds that race to truncate the same files.
    const int expected_ranks = argc == 3 && std::strcmp(argv[2], "1") == 0 ? 1 :
                               argc == 3 && std::strcmp(argv[2], "4") == 0 ? 4 : 0;
    if (!expected_ranks || ranks != expected_ranks) {
        if (rank == 0) {
            if (!expected_ranks)
                std::cerr << "Usage: " << argv[0] << " OUTPUT_PATH EXPECTED_RANKS(1|4)\n";
            else
                std::cerr << "MPI world size mismatch: expected " << expected_ranks
                          << ", got " << ranks << "; check the MPI compiler and launcher\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (rank == 0) std::cout << "DRM writer: verified MPI world size " << ranks << '\n';
    try {
        // Four ranks exercise a rank with zero points and interleaved views.
        const int active = ranks == 1 ? 1 : ranks - 1;
        std::vector<int> ids;
        for (int id = 0; id < 7; ++id)
            if (id % active == rank) ids.push_back(id);
        // Opening twice verifies truncation of an older, longer file.
        for (int frames : {9, 5}) {
            flat_output::FlatWriter writer;
            writer.open(MPI_COMM_WORLD, argv[1], 7, ids, 2);
            for (int frame = 0; frame < frames; ++frame) {
                double* out = writer.next_frame();
                for (size_t p = 0; p < ids.size(); ++p)
                    for (int axis = 0; axis < 3; ++axis)
                        out[p * 3 + axis] = frame * 1000 + ids[p] * 10 + axis + 0.25;
                writer.finish_frame();
            }
            writer.close(); // Flush one trailing frame.
            if (rank == 0) {
                std::ifstream input(argv[1], std::ios::binary | std::ios::ate);
                assert(input.tellg() == frames * 7 * 3 * static_cast<int>(sizeof(double)));
                input.seekg(0);
                for (int frame = 0; frame < frames; ++frame)
                    for (int id = 0; id < 7; ++id)
                        for (int axis = 0; axis < 3; ++axis) {
                            double value;
                            input.read(reinterpret_cast<char*>(&value), sizeof(value));
                            assert(value == frame * 1000 + id * 10 + axis + 0.25);
                        }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        // Simultaneous differently sized/rotated boxes, each with its own
        // sample rate and partial final batch. The payload is an affine field.
        drm_output::BoxGeometry boxes[2] = {
            drm_output::box_faces(drm_output::BoxSpec{{10,20,4,6,3,1,0}}),
            drm_output::box_faces(drm_output::BoxSpec{{30,50,40,60,30,1,37}})};
        flat_output::FlatWriter writers[2];
        std::vector<int> box_ids[2];
        const int rates[2] = {2, 3};
        for (int box = 0; box < 2; ++box) {
            for (auto& point : boxes[box].points)
                point[2] += 25.125 + 0.0125 * point[0] - 0.0025 * point[1];
            for (size_t id = 0; id < boxes[box].points.size(); ++id)
                if (static_cast<int>(id) % active == rank) box_ids[box].push_back(id);
            const std::string path = std::string(argv[1]) + ".box" + std::to_string(box);
            std::vector<drm_output::Material> local_materials;
            for (int id : box_ids[box]) {
                const auto& point = boxes[box].points[id];
                local_materials.push_back(drm_output::Material{{300 + point[2], 700 + point[2],
                    1800 + point[2], id == 0 ? NAN : id == 1 ? INFINITY : 25.125,
                    id == 1 ? INFINITY : 50.25}});
            }
            const auto materials = drm_output::gather_materials(MPI_COMM_WORLD,
                boxes[box].points.size(), box_ids[box], local_materials);
            if (rank == 0) drm_output::write_coordinates(path + ".coords", boxes[box].points, materials);
            else assert(materials.empty());
            writers[box].open(MPI_COMM_WORLD, path.c_str(), boxes[box].points.size(), box_ids[box], 2);
        }
        for (int step = 0; step < 14; ++step)
            for (int box = 0; box < 2; ++box) {
                if (step % rates[box]) continue;
                double* out = writers[box].next_frame();
                for (size_t p = 0; p < box_ids[box].size(); ++p)
                    for (int axis = 0; axis < 3; ++axis)
                        out[p*3+axis] = 10000*box + 1000*step + boxes[box].points[box_ids[box][p]][axis];
                writers[box].finish_frame();
            }
        for (int box = 0; box < 2; ++box) {
            writers[box].close();
            if (rank == 0) {
                const std::string path = std::string(argv[1]) + ".box" + std::to_string(box);
                std::ifstream input(path, std::ios::binary | std::ios::ate);
                const size_t frames = 13 / rates[box] + 1;
                assert(input.tellg() == static_cast<std::streamoff>(frames * boxes[box].points.size() * 3 * sizeof(double)));
                input.seekg(0);
                std::ifstream coordinates(path + ".coords");
                std::vector<drm_output::Point> points;
                std::string line;
                while (std::getline(coordinates, line)) {
                    if (!line.empty() && line[0] == '#') continue;
                    std::istringstream row(line);
                    drm_output::Point point;
                    std::string extra;
                    assert(row >> point[0] >> point[1] >> point[2]);
                    const drm_output::Material expected{{300 + point[2], 700 + point[2],
                        1800 + point[2], points.empty() ? NAN : points.size() == 1 ? INFINITY : 25.125,
                        points.size() == 1 ? INFINITY : 50.25}};
                    for (double value : expected) {
                        assert(row >> extra);
                        const double actual = std::stod(extra);
                        assert(std::isnan(value) ? std::isnan(actual) : actual == value);
                    }
                    assert(!(row >> extra));
                    points.push_back(point);
                }
                assert(points == boxes[box].points);
                for (int step = 0; step < 14; step += rates[box])
                    for (const auto& point : points)
                        for (int axis = 0; axis < 3; ++axis) {
                            double value;
                            input.read(reinterpret_cast<char*>(&value), sizeof(value));
                            assert(value == 10000*box + 1000*step + point[axis]);
                        }
            }
        }
        if (rank == 0) std::cout << "DRM MPI layout, empty-rank, batch, truncation and multiple-box regressions passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
