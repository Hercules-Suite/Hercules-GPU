/* Linux NUMA placement for immutable material tables. Setup only; no libnuma
 * dependency and no changes to CPU binding or the solver's later heap policy. */
#ifndef HERCULES_MATERIAL_NUMA_H
#define HERCULES_MATERIAL_NUMA_H

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#ifdef __linux__
#include <dirent.h>
#include <sched.h>
#include <linux/mempolicy.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

static inline bool material_numa_error(char *error, size_t capacity,
                                       const char *operation) {
    snprintf(error, capacity, "%s: %s", operation, strerror(errno));
    return false;
}

/* Every CPU allowed to this rank must belong to the same NUMA domain. Do not
 * infer locality from the CPU it happens to be executing on at one instant.
 * Dynamic masks handle sparse CPU numbering and hosts with >1024 CPUs. */
static inline bool material_numa_domain(int *node, char *error, size_t capacity,
                                        std::vector<int> *cpus = NULL) {
    if (cpus != NULL) cpus->clear();
#ifdef __linux__
    *node = -1;
    for (int ncpu = 128; ncpu <= (1 << 20); ncpu *= 2) {
        cpu_set_t *mask = CPU_ALLOC(ncpu);
        if (mask == NULL) return material_numa_error(error, capacity, "CPU_ALLOC");
        const size_t bytes = CPU_ALLOC_SIZE(ncpu);
        CPU_ZERO_S(bytes, mask);
        if (sched_getaffinity(0, bytes, mask) != 0) {
            const int saved_errno = errno;
            CPU_FREE(mask);
            if (saved_errno == EINVAL) continue;
            errno = saved_errno;
            return material_numa_error(error, capacity, "sched_getaffinity");
        }
        for (int cpu = 0; cpu < ncpu; ++cpu) {
            if (!CPU_ISSET_S(cpu, bytes, mask)) continue;
            char path[128];
            snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d", cpu);
            DIR *dir = opendir(path);
            if (dir == NULL) {
                const int saved_errno = errno;
                CPU_FREE(mask);
                errno = saved_errno;
                return material_numa_error(error, capacity, "open CPU NUMA topology");
            }
            int cpu_node = -1;
            struct dirent *entry;
            while ((entry = readdir(dir)) != NULL) {
                if (strncmp(entry->d_name, "node", 4) != 0) continue;
                char *end = NULL;
                errno = 0;
                const long id = strtol(entry->d_name + 4, &end, 10);
                if (errno == 0 && end != entry->d_name + 4 && *end == '\0' &&
                    id >= 0 && id < (1 << 20)) cpu_node = (int)id;
            }
            closedir(dir);
            if (cpu_node < 0 || (*node >= 0 && cpu_node != *node)) {
                CPU_FREE(mask);
                snprintf(error, capacity, cpu_node < 0
                    ? "CPU %d has no readable Linux NUMA topology"
                    : "CPU affinity spans NUMA domains; bind each rank (and its threads) to one domain; conflicting CPU %d",
                    cpu);
                return false;
            }
            *node = cpu_node;
            if (cpus != NULL) cpus->push_back(cpu);
        }
        CPU_FREE(mask);
        if (*node >= 0) return true;
        break;
    }
    snprintf(error, capacity, "empty or unsupported CPU affinity mask");
#else
    (void)node;
    snprintf(error, capacity, "NUMA material storage requires Linux CPU affinity, sysfs and memory-policy syscalls");
#endif
    return false;
}

struct material_numa_policy_t {
    int saved_mode;
    std::vector<unsigned long> saved_mask;
    bool active;
    material_numa_policy_t() : saved_mode(0), active(false) {}

    bool begin(int node, char *error, size_t capacity) {
#ifdef __linux__
        const unsigned long word_bits = 8 * sizeof(unsigned long);
        for (unsigned long bits = 128; bits <= (1UL << 20); bits *= 2) {
            if (node < 0 || (unsigned long)node >= bits) continue;
            saved_mask.assign((bits + word_bits - 1) / word_bits, 0);
            if (syscall(SYS_get_mempolicy, &saved_mode, saved_mask.data(),
                        bits, NULL, 0UL) != 0) {
                if (errno == EINVAL) continue;
                return material_numa_error(error, capacity, "get_mempolicy");
            }
            std::vector<unsigned long> target(saved_mask.size(), 0);
            target[node / word_bits] |= 1UL << (node % word_bits);
            if (syscall(SYS_set_mempolicy, MPOL_BIND | MPOL_F_STATIC_NODES,
                        target.data(), bits) != 0)
                return material_numa_error(error, capacity, "set_mempolicy(MPOL_BIND)");
            active = true;
            return true;
        }
        snprintf(error, capacity, "unsupported NUMA node mask size");
#else
        (void)node;
        snprintf(error, capacity, "NUMA memory policy requires Linux");
#endif
        return false;
    }

    bool restore(char *error, size_t capacity) {
        if (!active) return true;
#ifdef __linux__
        const unsigned long bits = saved_mask.size() * 8 * sizeof(unsigned long);
        if (syscall(SYS_set_mempolicy, saved_mode,
                    saved_mode == MPOL_DEFAULT ? NULL : saved_mask.data(),
                    saved_mode == MPOL_DEFAULT ? 0UL : bits) != 0)
            return material_numa_error(error, capacity, "restore memory policy");
#else
        snprintf(error, capacity, "NUMA memory policy requires Linux");
        return false;
#endif
        active = false;
        return true;
    }
};

/* Fault each page on its NUMA-local leader before MPI receives can touch it.
 * Stores are overwritten by the complete model transfer before publication.
 * Checking at most 32 addresses per array detects incompatible backing or
 * placement policies without scanning every page with a syscall. */
static inline bool material_numa_touch(void *base, size_t bytes, int node,
                                      int *checked, char *error, size_t capacity) {
#ifdef __linux__
    const long page_size = sysconf(_SC_PAGESIZE);
    if (page_size <= 0 || bytes == 0) {
        snprintf(error, capacity, "invalid page or allocation size");
        return false;
    }
    volatile unsigned char *data = (volatile unsigned char *)base;
    for (size_t offset = 0; offset < bytes; offset += (size_t)page_size)
        data[offset] = 0;
    data[bytes - 1] = 0;
    const size_t samples = ((uintptr_t)base % (size_t)page_size + bytes - 1)
                         / (size_t)page_size + 1;
    const size_t count = samples < 32 ? samples : 32;
    for (size_t i = 0; i < count; ++i) {
        const size_t offset = count == 1 ? 0 : i == count - 1 ? bytes - 1 :
            ((bytes - 1) / (count - 1)) * i;
        int actual = -1;
        if (syscall(SYS_get_mempolicy, &actual, NULL, 0UL,
                    (unsigned char *)base + offset,
                    (unsigned long)(MPOL_F_NODE | MPOL_F_ADDR)) != 0)
            return material_numa_error(error, capacity, "query material page NUMA node");
        if (actual != node) {
            snprintf(error, capacity, "material page is on NUMA node %d, expected %d; check MPI shared-memory backing and job memory policy",
                     actual, node);
            return false;
        }
        ++*checked;
    }
    return true;
#else
    (void)base; (void)bytes; (void)node; (void)checked;
    snprintf(error, capacity, "NUMA page placement requires Linux");
    return false;
#endif
}
#endif
