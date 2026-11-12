/* Setup-only survivor and CPU placement for a NUMA material-model handoff.
 * No MPI/CUDA dependency: the caller supplies observed rank/GPU topology. */
#ifndef HERCULES_UNIFIED_NUMA_H
#define HERCULES_UNIFIED_NUMA_H
#include "material_numa.h"
#include <algorithm>
#include <map>

struct unified_numa_plan_t {
    std::vector<int> ranks; /* Original node-local rank for each target GPU ordinal. */
    std::vector<int> nodes;
    std::vector<std::vector<int> > cpus;
};

inline bool unified_numa_plan(const std::vector<int>& rank_nodes,
        const std::vector<std::vector<int> >& rank_cpus,
        const std::vector<int>& gpu_nodes, unified_numa_plan_t *plan,
        char *error, size_t capacity) {
    struct domain_t {
        std::vector<int> ranks, cpus, slots;
    };
    std::map<int, domain_t> domains;
    std::map<int, int> cpu_nodes;
    if (rank_nodes.empty() || rank_cpus.size() != rank_nodes.size() ||
        gpu_nodes.empty() || gpu_nodes.size() > rank_nodes.size()) {
        snprintf(error, capacity, "invalid NUMA handoff rank/GPU counts");
        return false;
    }
    for (size_t r = 0; r < rank_nodes.size(); ++r) {
        if (rank_nodes[r] < 0 || rank_cpus[r].empty()) {
            snprintf(error, capacity, "missing NUMA affinity for local rank %zu", r);
            return false;
        }
        domain_t& domain = domains[rank_nodes[r]];
        domain.ranks.push_back((int)r);
        for (int cpu : rank_cpus[r]) {
            if (cpu < 0 || cpu >= (1 << 20) ||
                (cpu_nodes.count(cpu) && cpu_nodes[cpu] != rank_nodes[r])) {
                snprintf(error, capacity, "inconsistent NUMA affinity for CPU %d", cpu);
                return false;
            }
            cpu_nodes[cpu] = rank_nodes[r];
            domain.cpus.push_back(cpu);
        }
    }
    for (auto& item : domains) {
        std::vector<int>& cpus = item.second.cpus;
        std::sort(cpus.begin(), cpus.end());
        cpus.erase(std::unique(cpus.begin(), cpus.end()), cpus.end());
    }
    plan->ranks.assign(gpu_nodes.size(), -1);
    plan->nodes.assign(gpu_nodes.size(), -1);
    plan->cpus.assign(gpu_nodes.size(), std::vector<int>());
    /* Reserve known GPU domains first, then balance unknown-GPU counts over
     * available domains (ties use the first original rank). Group their GPU
     * ordinals into contiguous domain blocks below. This preserves material
     * locality but cannot establish GPU locality without reported topology. */
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t slot = 0; slot < gpu_nodes.size(); ++slot) {
            const int gpu_node = gpu_nodes[slot];
            if ((gpu_node >= 0) != (pass == 0)) continue;
            auto chosen = domains.end();
            if (gpu_node >= 0) {
                chosen = domains.find(gpu_node);
            } else {
                for (auto it = domains.begin(); it != domains.end(); ++it) {
                    const domain_t& d = it->second;
                    if (d.slots.size() >= d.ranks.size() ||
                        d.slots.size() >= d.cpus.size()) continue;
                    if (chosen == domains.end() ||
                        d.slots.size() < chosen->second.slots.size() ||
                        (d.slots.size() == chosen->second.slots.size() &&
                         d.ranks.front() < chosen->second.ranks.front())) chosen = it;
                }
            }
            if (chosen == domains.end() ||
                chosen->second.slots.size() >= chosen->second.ranks.size() ||
                chosen->second.slots.size() >= chosen->second.cpus.size()) {
                snprintf(error, capacity,
                    "GPU %zu (NUMA node %d) has no available rank/CPU in its domain; cover GPU domains with Phase-A bindings",
                    slot, gpu_node);
                return false;
            }
            domain_t& d = chosen->second;
            plan->ranks[slot] = d.ranks[d.slots.size()];
            plan->nodes[slot] = chosen->first;
            d.slots.push_back((int)slot);
        }
    }
    /* Keep the balanced, capacity-limited domain counts, but avoid alternating
     * GPU ordinals between domains. On a two-domain/four-GPU node, unknown
     * GPUs 0/1 belong to the first domain and 2/3 to the second. Known GPU
     * assignments remain authoritative, even when they interrupt these blocks. */
    std::vector<size_t> unknown_slots;
    std::vector<int> unknown_ranks;
    for (size_t slot = 0; slot < gpu_nodes.size(); ++slot) {
        if (gpu_nodes[slot] < 0) {
            unknown_slots.push_back(slot);
            unknown_ranks.push_back(plan->ranks[slot]);
        }
    }
    std::sort(unknown_ranks.begin(), unknown_ranks.end(), [&](int a, int b) {
        const int first_a = domains.at(rank_nodes[a]).ranks.front();
        const int first_b = domains.at(rank_nodes[b]).ranks.front();
        return first_a != first_b ? first_a < first_b : a < b;
    });
    for (size_t i = 0; i < unknown_slots.size(); ++i) {
        plan->ranks[unknown_slots[i]] = unknown_ranks[i];
        plan->nodes[unknown_slots[i]] = rank_nodes[unknown_ranks[i]];
    }
    for (auto& item : domains) item.second.slots.clear();
    for (size_t slot = 0; slot < plan->nodes.size(); ++slot)
        domains[plan->nodes[slot]].slots.push_back((int)slot);
    for (auto& item : domains) {
        domain_t& d = item.second;
        std::sort(d.slots.begin(), d.slots.end());
        for (size_t i = 0; i < d.slots.size(); ++i) {
            const size_t first = i * d.cpus.size() / d.slots.size();
            const size_t last = (i + 1) * d.cpus.size() / d.slots.size();
            plan->cpus[d.slots[i]].assign(d.cpus.begin() + first,
                                        d.cpus.begin() + last);
        }
    }
    return true;
}

/* CUDA ordinals honor CUDA_VISIBLE_DEVICES. Resolve the selected ordinal's
 * PCI bus ID to a Linux NUMA domain; -1 means unavailable, not domain zero. */
inline int unified_gpu_numa_node(const char *bus_id) {
#ifdef __linux__
    unsigned domain, bus, device, function;
    char extra;
    if (sscanf(bus_id, "%x:%x:%x.%x%c", &domain, &bus, &device, &function,
               &extra) != 4 || domain > 0xffff || bus > 0xff ||
        device > 0x1f || function > 7) return -1;
    char path[128];
    snprintf(path, sizeof(path),
        "/sys/bus/pci/devices/%04x:%02x:%02x.%x/numa_node",
        domain, bus, device, function);
    FILE *file = fopen(path, "r");
    if (file == NULL) return -1;
    int node = -1;
    if (fscanf(file, "%d", &node) != 1 || node < 0 || node >= (1 << 20)) node = -1;
    fclose(file);
    return node;
#else
    (void)bus_id;
    return -1;
#endif
}

inline bool unified_numa_bind(const std::vector<int>& cpus,
                                     char *error, size_t capacity) {
#ifdef __linux__
    if (cpus.empty() || cpus.back() < 0 || cpus.back() >= (1 << 20)) {
        snprintf(error, capacity, "invalid Phase-B NUMA CPU mask");
        return false;
    }
    const size_t bytes = CPU_ALLOC_SIZE(cpus.back() + 1);
    cpu_set_t *mask = CPU_ALLOC(cpus.back() + 1);
    if (mask == NULL) return material_numa_error(error, capacity, "CPU_ALLOC");
    CPU_ZERO_S(bytes, mask);
    for (int cpu : cpus) {
        if (cpu < 0 || cpu > cpus.back()) {
            CPU_FREE(mask);
            snprintf(error, capacity, "invalid Phase-B NUMA CPU index");
            return false;
        }
        CPU_SET_S(cpu, bytes, mask);
    }
    const int result = sched_setaffinity(0, bytes, mask);
    const int saved_errno = errno;
    CPU_FREE(mask);
    if (result == 0) {
        int node;
        std::vector<int> actual;
        if (!material_numa_domain(&node, error, capacity, &actual)) return false;
        if (actual == cpus) return true;
        snprintf(error, capacity,
            "Phase-B CPU mask was restricted by the job cpuset; allow each survivor its domain CPU share");
        return false;
    }
    errno = saved_errno;
    return material_numa_error(error, capacity, "bind Phase-B NUMA CPU mask");
#else
    (void)cpus;
    snprintf(error, capacity, "NUMA CPU binding requires Linux");
    return false;
#endif
}
#endif
