/* Native constitutive checks for the two RIVA-Sand 02 input interfaces.
 * Parser/storage/dispatch coverage lives in docs-site's source-wiring tests;
 * run_suite.sh also requires full CPU and CUDA reference-row CSV equivalence.
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

static void require(bool ok, const char *message)
{
    if (!ok) throw std::runtime_error(message);
}

static bool same(const riva_ib_state_t &a, const riva_ib_state_t &b)
{
    double av[RIVA_IB_STATE_VALUE_COUNT], bv[RIVA_IB_STATE_VALUE_COUNT];
    require(riva_ib_state_values(&a, av) == RIVA_IB_STATE_VALUE_COUNT &&
            riva_ib_state_values(&b, bv) == RIVA_IB_STATE_VALUE_COUNT,
            "wrong state width");
    for (int i = 0; i < RIVA_IB_STATE_VALUE_COUNT; ++i) {
        require(std::isfinite(av[i]) && std::isfinite(bv[i]), "nonfinite state");
        if (av[i] != bv[i]) return false;
    }
    return a.base.initialized == b.base.initialized &&
           a.base.geostatic_admitted == b.base.geostatic_admitted;
}

static riva_ib_state_t initialize(const riva_ib_parameters_t &p,
                                 const riva_material_parameters_t &m)
{
    const tensor_t stress = {-19.4, -19.4, -40.0, 10.0, 0.0, 0.0};
    tensor_t reference = stress;
    reference.xy = 0.0;
    riva_ib_state_t s = {};
    const double e = riva_void_ratio_from_material_relative_density(
        &m, 0.66296296296296309);
    require(riva_material_parameters_valid(&p.base, &m), "invalid custom row");
    require(riva_ib_initialize_material(&p, &m, stress, e, &s) &&
            riva_ib_begin_dynamic_phase(&p, &m, &reference, &s),
            "custom initialization failed");
    require(s.base.void_ratio == e, "custom void-ratio bounds ignored");
    return s;
}

static void advance(const riva_ib_parameters_t &p,
                    const riva_material_parameters_t &m,
                    riva_ib_state_t &s, int step, int nsub)
{
    const double angle = 2.0*3.14159265358979323846/32.0;
    const double shear = 0.002*(std::sin((step+1)*angle)-std::sin(step*angle));
    const tensor_t deps = {0.0, 0.0, 0.0, shear, 0.0, 0.0};
    tensor_t stress = s.base.stress;
    riva_update_info_t info = {};
    require(riva_ib_update_material(&p, &m, deps, nsub, &s, &stress, &info),
            "custom update failed");
}

int main()
{
    try {
        require(riva_ib_checkpoint_revision(0) == 3u &&
                riva_ib_checkpoint_revision(1) == 4u,
                "custom density correction must reject old Custom02 checkpoints");
        const riva_ib_parameters_t p = riva_ib_reference_parameters(1.0);
        const riva_material_parameters_t frozen =
            riva_reference_material_parameters(&p.base);
        riva_material_parameters_t reference = {};
        reference.M=1.25; reference.kd=1.125;
        reference.h=122.44207260468994; reference.m=0.945;
        reference.zeta=0.025; reference.e_max=0.78; reference.e_min=0.51;
        reference.Q=10.0; reference.R=1.5; reference.n_G=0.65;
        const std::array<double riva_material_parameters_t::*, 10> fields = {
            &riva_material_parameters_t::M, &riva_material_parameters_t::kd,
            &riva_material_parameters_t::h, &riva_material_parameters_t::m,
            &riva_material_parameters_t::zeta, &riva_material_parameters_t::e_max,
            &riva_material_parameters_t::e_min, &riva_material_parameters_t::Q,
            &riva_material_parameters_t::R, &riva_material_parameters_t::n_G};
        for (auto field : fields)
            require(reference.*field == frozen.*field, "reference value drift");
        for (const int nsub : {1, 4, 16}) {
            for (auto field : fields) {
                auto custom = reference;
                custom.*field *= 1.05;
                auto a = initialize(p, frozen), b = initialize(p, reference);
                auto c = initialize(p, custom), restarted = c;
                bool changed = !same(a, c);
                for (int step = 0; step < 160; ++step) {
                    advance(p, custom, c, step, nsub);
                    advance(p, frozen, a, step, nsub);
                    advance(p, reference, b, step, nsub);
                    advance(p, custom, restarted, step, nsub);
                    require(same(a, b), "default/custom reference mismatch");
                    require(same(c, restarted), "custom point/restart mismatch");
                    changed |= !same(a, c);
                    if (step == 48) {
                        unsigned char checkpoint[sizeof(c)];
                        std::memcpy(checkpoint, &c, sizeof(c));
                        restarted = {};
                        std::memcpy(&restarted, checkpoint, sizeof(c));
                    }
                }
                require(changed, "custom physical input did not affect history");
            }
        }
        auto invalid = reference;
        invalid.e_min = 0.7;
        require(!riva_material_parameters_valid(&p.base, &invalid),
                "out-of-range shakedown anchor accepted");
        std::cout << "PASS: RIVASand02/RIVASandCustom02 exact reference equality, "
                     "all ten physical inputs active, per-material isolation, "
                     "custom byte restart, and bounds (1/4/16 substeps)\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
