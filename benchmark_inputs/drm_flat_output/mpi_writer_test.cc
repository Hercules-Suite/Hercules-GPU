// Real collective I/O regression; no CUDA or solver build required.
#include "../../quake/forward/drm_output_mpi.h"
#include <cassert>
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (argc != 2) MPI_Abort(MPI_COMM_WORLD, 1);
    try {
        // Four ranks exercise a rank with zero points and interleaved views.
        const int active = ranks == 1 ? 1 : ranks - 1;
        std::vector<int> ids;
        for (int id = 0; id < 7; ++id)
            if (id % active == rank) ids.push_back(id);
        // Opening twice verifies truncation of an older, longer file.
        for (int frames : {9, 5}) {
            drm_output::FlatWriter writer;
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
        if (rank == 0) std::cout << "DRM MPI layout, empty-rank, batch, truncation regressions passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
