/* Rank-independent flat [snapshot][point][XYZ] native float64 output. */
#ifndef DRM_OUTPUT_MPI_H
#define DRM_OUTPUT_MPI_H
#include <mpi.h>
#include <climits>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace drm_output {
inline void mpi_check(int code, const char* operation) {
    if (code == MPI_SUCCESS) return;
    char error[MPI_MAX_ERROR_STRING];
    int length = 0;
    MPI_Error_string(code, error, &length);
    throw std::runtime_error(std::string(operation) + ": " +
                             std::string(error, length));
}

class FlatWriter {
    MPI_File file_ = MPI_FILE_NULL;
    MPI_Datatype view_ = MPI_DATATYPE_NULL;
    MPI_Offset snapshot_bytes_ = 0;
    MPI_Offset frames_ = 0;
    int values_ = 0, capacity_ = 0, buffered_ = 0;
    std::vector<double> buffer_;
public:
    void open(MPI_Comm comm, const char* path, int total_points,
              const std::vector<int>& ids, int capacity) {
        if (total_points <= 0 || capacity <= 0 ||
            ids.size() > static_cast<size_t>(INT_MAX / 3 / capacity))
            throw std::runtime_error("DRM MPI batch exceeds int-count limits; reduce drm_output_batch_size");
        const int64_t snapshot_bytes = static_cast<int64_t>(total_points) * 3 * sizeof(double);
        if (snapshot_bytes > std::numeric_limits<MPI_Aint>::max() ||
            snapshot_bytes > std::numeric_limits<MPI_Offset>::max())
            throw std::runtime_error("DRM snapshot exceeds MPI address/offset limits");
        snapshot_bytes_ = static_cast<MPI_Offset>(snapshot_bytes);
        values_ = static_cast<int>(ids.size()) * 3;
        capacity_ = capacity;
        buffer_.resize(values_ ? static_cast<size_t>(values_) * capacity : 1);
        if (!ids.empty()) {
            std::vector<int> lengths(ids.size(), 3);
            std::vector<MPI_Aint> offsets(ids.size());
            for (size_t i = 0; i < ids.size(); ++i) {
                if (ids[i] < 0 || ids[i] >= total_points || (i && ids[i] <= ids[i-1]))
                    throw std::runtime_error("DRM point IDs must be increasing and in range");
                offsets[i] = static_cast<MPI_Aint>(ids[i]) * 3 * sizeof(double);
            }
            MPI_Datatype raw;
            mpi_check(MPI_Type_create_hindexed(static_cast<int>(ids.size()),
                       lengths.data(), offsets.data(), MPI_DOUBLE, &raw), "DRM file layout");
            mpi_check(MPI_Type_create_resized(raw, 0, snapshot_bytes_, &view_), "DRM frame extent");
            mpi_check(MPI_Type_free(&raw), "DRM raw type free");
        } else {
            // Empty ranks participate in every collective with count zero.
            mpi_check(MPI_Type_contiguous(1, MPI_DOUBLE, &view_), "DRM empty view");
        }
        mpi_check(MPI_Type_commit(&view_), "DRM view commit");
        mpi_check(MPI_File_open(comm, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                               MPI_INFO_NULL, &file_), "DRM file open");
        mpi_check(MPI_File_set_errhandler(file_, MPI_ERRORS_RETURN), "DRM error handler");
        mpi_check(MPI_File_set_size(file_, 0), "DRM file truncate");
        mpi_check(MPI_File_set_view(file_, 0, MPI_DOUBLE, view_, "native", MPI_INFO_NULL),
                  "DRM file view");
    }
    double* next_frame() { return buffer_.data() + static_cast<size_t>(buffered_) * values_; }
    void finish_frame() { if (++buffered_ == capacity_) flush(); }
    void flush() {
        if (!buffered_) return;
        if (frames_ > std::numeric_limits<MPI_Offset>::max() / snapshot_bytes_ - buffered_)
            throw std::runtime_error("DRM file exceeds MPI offset limits");
        MPI_Status status;
        // Explicit offsets are etypes within this rank's view, not bytes.
        mpi_check(MPI_File_write_at_all(file_, frames_ * values_, buffer_.data(),
                                       buffered_ * values_, MPI_DOUBLE, &status), "DRM write");
        int written = 0;
        mpi_check(MPI_Get_count(&status, MPI_DOUBLE, &written), "DRM write count");
        if (written != buffered_ * values_) throw std::runtime_error("Short DRM write");
        frames_ += buffered_;
        buffered_ = 0;
    }
    void close() {
        if (file_ == MPI_FILE_NULL) return;
        flush();
        mpi_check(MPI_File_set_size(file_, frames_ * snapshot_bytes_), "DRM final size");
        mpi_check(MPI_File_close(&file_), "DRM close");
        mpi_check(MPI_Type_free(&view_), "DRM view free");
        std::vector<double>().swap(buffer_);
    }
};
} // namespace drm_output
#endif
