#ifndef HERCULES_CPML_PARTITION_MODEL_H
#define HERCULES_CPML_PARTITION_MODEL_H

/* Setup/finalization-only cost model. No CUDA, MPI, or physical state here. */
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>

namespace cpml_partition {

static constexpr double max_auto_weight = 16.0;
static constexpr double max_relative_error = 0.10;
static constexpr double min_fraction_span = 0.05;
static constexpr double max_numeric_weight =
    static_cast<double>(std::numeric_limits<int32_t>::max()) / 100.0;

/* Unlike the historical optional-string reader, accept indented keys too.
 * These controls take exactly one whitespace-free token. */
inline int read_control(FILE* file, const char* key, std::string& value) {
  std::rewind(file);
  char* buffer = nullptr;
  size_t capacity = 0;
  int found = 0;
  while (getline(&buffer, &capacity, file) >= 0) {
    std::string line(buffer);
    line = line.substr(0, line.find('#'));
    const size_t equal = line.find('=');
    if (equal == std::string::npos) continue;
    std::istringstream left(line.substr(0, equal)), right(line.substr(equal + 1));
    std::string name, extra, token;
    left >> name;
    if (name != key || (left >> extra)) continue;
    if (found || !(right >> token) || (right >> extra)) {
      found = -1;
      break;
    }
    value = token;
    found = 1;
  }
  if (std::ferror(file)) found = -1;
  std::free(buffer);
  return found;
}

inline bool parse_weight(const std::string& text, double& weight, bool& automatic) {
  std::string lower = text;
  for (char& c : lower)
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
  automatic = lower == "auto";
  if (automatic) {
    weight = 1.0;
    return true;
  }
  char* end = nullptr;
  errno = 0;
  const double value = std::strtod(text.c_str(), &end);
  if (errno || end == text.c_str() || *end != '\0' ||
      !std::isfinite(value) || value < 1.0 || value > max_numeric_weight)
    return false;
  weight = value;
  return true;
}

struct Sample {
  double elements;
  double pml;
  double seconds;
};

struct Estimate {
  bool valid = false;
  double weight = 1.0;
  double ordinary_cost = 0.0;
  double pml_cost = 0.0;
  double relative_error = 0.0;
  size_t samples = 0;
  const char* reason = "insufficient samples";
};

/* T/N = ordinary_cost + (pml_cost - ordinary_cost) * P/N. Using all
 * ranks, rather than unrelated timer extrema, also handles unequal N after
 * weighting. A scalar fit can be unsuitable for mixed physics, nonlinear
 * substepping, or heterogeneous material work: reject poorly explained data.
 * Costs include existing Compute Physics timing, not source I/O or standalone
 * MPI exchange timers. They are effective costs, not isolated kernel timings. */
inline Estimate estimate(const std::vector<Sample>& samples) {
  Estimate out;
  double mx = 0.0, my = 0.0, xmin = 1.0, xmax = 0.0;
  for (const Sample& s : samples) {
    if (!std::isfinite(s.elements) || !std::isfinite(s.pml) ||
        !std::isfinite(s.seconds) || s.elements < 0.0 || s.pml < 0.0 ||
        s.pml > s.elements || s.seconds < 0.0) {
      out.reason = "invalid counts or timings";
      return out;
    }
    if (s.elements == 0.0) continue;
    if (s.seconds < 0.1) {
      out.reason = "compute interval too short";
      return out;
    }
    const double x = s.pml / s.elements, y = s.seconds / s.elements;
    ++out.samples;
    mx += (x - mx) / out.samples;
    my += (y - my) / out.samples;
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
  }
  if (out.samples < 4) return out;
  if (xmax - xmin < min_fraction_span) {
    out.reason = "insufficient variation in PML fractions";
    return out;
  }
  double xx = 0.0, xy = 0.0;
  for (const Sample& s : samples) {
    if (s.elements == 0.0) continue;
    const double dx = s.pml / s.elements - mx;
    xx += dx * dx;
    xy += dx * (s.seconds / s.elements - my);
  }
  const double slope = xy / xx;
  out.ordinary_cost = my - slope * mx;
  out.pml_cost = out.ordinary_cost + slope;
  if (!(out.ordinary_cost > 0.0) || !(out.pml_cost > 0.0) ||
      !std::isfinite(out.ordinary_cost) || !std::isfinite(out.pml_cost)) {
    out.reason = "nonpositive or nonfinite fitted costs";
    return out;
  }
  double residual = 0.0, energy = 0.0;
  for (const Sample& s : samples) {
    if (s.elements == 0.0) continue;
    const double predicted = out.ordinary_cost * (s.elements - s.pml) +
                             out.pml_cost * s.pml;
    residual += (predicted - s.seconds) * (predicted - s.seconds);
    energy += s.seconds * s.seconds;
  }
  out.relative_error = std::sqrt(residual / energy);
  if (!std::isfinite(out.relative_error) || out.relative_error > max_relative_error) {
    out.reason = "scalar cost model does not explain rank timings";
    return out;
  }
  out.weight = std::max(1.0, std::min(max_auto_weight,
                                    out.pml_cost / out.ordinary_cost));
  out.valid = true;
  out.reason = "accepted";
  return out;
}

/* Stable across host byte orders. Used only as a conservative cache identity,
 * not as a cryptographic digest or mesh-equivalence check. */
inline uint64_t hash(const std::string& text, uint64_t value = UINT64_C(14695981039346656037)) {
  for (unsigned char c : text) {
    value ^= c;
    value *= UINT64_C(1099511628211);
  }
  return value;
}

inline bool load(const std::string& path, uint64_t key, double& weight) {
  std::ifstream file(path);
  std::string magic, trailing;
  uint64_t stored_key = 0;
  double stored_weight = 0.0, error = 0.0;
  size_t count = 0;
  if (!(file >> magic >> std::hex >> stored_key >> std::dec >> stored_weight >> count >> error) ||
      magic != "HERCULES_CPML_COST_V1" || stored_key != key ||
      !std::isfinite(stored_weight) || stored_weight < 1.0 ||
      stored_weight > max_auto_weight || count < 4 || !std::isfinite(error) ||
      error < 0.0 || error > max_relative_error || (file >> trailing))
    return false;
  weight = stored_weight;
  return true;
}

/* Atomic replacement: concurrent jobs may replace each other's complete
 * records, but readers never see a partially written calibration. A cache
 * holds one configuration; a mismatched key is a miss. Failure is nonfatal. */
inline bool save(const std::string& path, uint64_t key, const Estimate& fit) {
  if (!fit.valid || !std::isfinite(fit.weight) || fit.weight < 1.0 ||
      fit.weight > max_auto_weight || fit.samples < 4 ||
      !std::isfinite(fit.relative_error) || fit.relative_error < 0.0 ||
      fit.relative_error > max_relative_error)
    return false;
  // A mistaken cache path must not replace a mesh, input, or symlink.
  struct stat existing;
  if (lstat(path.c_str(), &existing) == 0) {
    if (!S_ISREG(existing.st_mode)) return false;
    std::ifstream old(path);
    std::string magic;
    if (!(old >> magic) || magic != "HERCULES_CPML_COST_V1") return false;
  } else if (errno != ENOENT) {
    return false;
  }
  std::string name = path + ".tmp.XXXXXX";
  std::vector<char> temp(name.begin(), name.end());
  temp.push_back('\0');
  const int fd = mkstemp(temp.data());
  if (fd < 0) return false;
  FILE* file = fdopen(fd, "w");
  if (!file) {
    close(fd);
    unlink(temp.data());
    return false;
  }
  bool ok = std::fprintf(file, "HERCULES_CPML_COST_V1 %016llx %.17g %zu %.17g\n",
      static_cast<unsigned long long>(key), fit.weight, fit.samples,
      fit.relative_error) > 0;
  if (std::fclose(file) != 0) ok = false;
  if (ok && std::rename(temp.data(), path.c_str()) != 0) ok = false;
  if (!ok) unlink(temp.data());
  return ok;
}

}  // namespace cpml_partition
#endif
