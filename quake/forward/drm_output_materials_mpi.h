/* Setup-only material gather in DRM sample order; no GPU field exchange. */
#ifndef DRM_OUTPUT_MATERIALS_MPI_H
#define DRM_OUTPUT_MATERIALS_MPI_H
#include "drm_output_coordinates.h"
#include "flat_output_mpi.h"

namespace drm_output {
// The caller has collectively verified exactly one owner per sample.
inline std::vector<Material> gather_materials(MPI_Comm comm, int total_points,
        const std::vector<int>& ids, const std::vector<Material>& local) {
    if (total_points <= 0 || ids.size() != local.size())
        throw std::runtime_error("Invalid DRM material gather dimensions");
    for (size_t i = 0; i < ids.size(); ++i)
        if (ids[i] < 0 || ids[i] >= total_points || (i && ids[i] <= ids[i-1]))
            throw std::runtime_error("DRM material IDs must be increasing and in range");
    int rank;
    flat_output::mpi_check(MPI_Comm_rank(comm, &rank), "DRM material rank");
    std::vector<Material> result(rank == 0 ? total_points : 0);
    static_assert(sizeof(Material) == 5 * sizeof(double), "Material rows must be packed");
    // Bound temporary storage and MPI counts even for very large boxes.
    const int chunk_points = 4096;
    std::vector<Material> buffer(chunk_points);
    size_t cursor = 0;
    for (int start = 0; start < total_points;) {
        const int count = std::min(chunk_points, total_points - start);
        std::fill(buffer.begin(), buffer.end(), Material{{0, 0, 0, 0, 0}});
        while (cursor < ids.size() && ids[cursor] < start + count) {
            buffer[ids[cursor] - start] = local[cursor];
            ++cursor;
        }
        // Nonowners contribute zero; the unique owner's NaN remains NaN.
        flat_output::mpi_check(MPI_Reduce(buffer.data(), rank == 0 ? result.data() + start : NULL,
            count * 5, MPI_DOUBLE, MPI_SUM, 0, comm), "DRM material gather");
        start += count;
    }
    return result;
}
} // namespace drm_output
#endif
