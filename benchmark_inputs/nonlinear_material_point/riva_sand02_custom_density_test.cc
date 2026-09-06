/* Custom bounds must preserve the e=0.601 shakedown anchor.
 * --probe emits material-specific diagnostics for independent Python checks.
 */
#define RIVA_STANDALONE
#include "../../quake/forward/nonlinear_riva_sand_intermediate_bias_research.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv)
{
    const bool probe = argc == 2 && std::string(argv[1]) == "--probe";
    try {
        const riva_ib_parameters_t p = riva_ib_reference_parameters(1.0);
        const double bounds[][2] = {
            {0.78, 0.51}, {0.819, 0.51}, {0.78, 0.5355}, {0.819, 0.5355},
            {0.7995, 0.52275}}; // also a midpoint of two interpolated rows
        if (probe)
            std::cout << "e_max,e_min,e0,reference_density,initial_density,raw,"
                         "shakedown,shear,density_weight,dense_bias_gate,"
                         "dense_transition,branch_compliance,phase_hardening\n";
        for (const auto &bound : bounds) {
            auto m = riva_reference_material_parameters(&p.base);
            m.e_max = bound[0]; m.e_min = bound[1];
            for (double e0 : {0.601, 0.581, 0.621}) {
                const tensor_t stress = {-19.4, -19.4, -40.0, 10.0, 0.0, 0.0};
                tensor_t reference = stress;
                reference.xy = 0.0;
                riva_ib_state_t s = {};
                if (!riva_ib_initialize_material(&p, &m, stress, e0, &s) ||
                    !riva_ib_begin_dynamic_phase(&p, &m, &reference, &s))
                    throw std::runtime_error("initialization failed");
                s.base.amplitude_reversals = 4;
                s.base.cyclic_amplitude = 0.6;
                s.base.cyclic_direction = {0, 0, 0, 1.0/std::sqrt(2.0), 0, 0};
                const double pressure = riva_pressure(s.base.stress);
                const double raw = riva_ib_raw_shakedown(&p, &m, &s, pressure);
                const double shake = riva_ib_shakedown(&p, &m, &s, pressure);
                const double shear = riva_ib_shakedown_shear(&p, &m, &s, pressure);
                const double weight = riva_ib_density_weight(&p, &m, &s);
                const double dense_gate = riva_ib_phase_dense_bias_gate(&p, &m, &s);
                const double transition = riva_ib_phase_dense_transition(&p, &m, &s);
                const double ref_density = riva_ib_reference_relative_density(&p, &m);
                if (probe) {
                    std::cout << std::setprecision(17)
                        << m.e_max << ',' << m.e_min << ',' << e0 << ','
                        << ref_density << ',' << riva_ib_initial_density(&s) << ','
                        << raw << ',' << shake << ',' << shear << ',' << weight
                        << ',' << dense_gate << ',' << transition << ','
                        << riva_ib_branch_compliance_multiplier(&p, &m, &s) << ','
                        << riva_ib_phase_bias_hardening(&p, &m, &s) << '\n';
                } else if (e0 == p.base.state_shakedown_reference_void_ratio) {
                    if (raw != 1.0 || shake != 1.0 || shear != 1.0 ||
                        weight != 0.0 || std::fabs(transition) > 1e-28) {
                        std::cerr << std::setprecision(17)
                            << "bounds=" << m.e_max << ',' << m.e_min
                            << " reference_density=" << ref_density
                            << " initial_density=" << riva_ib_initial_density(&s)
                            << " raw=" << raw << " shakedown=" << shake << '\n';
                        throw std::runtime_error("custom density moved shakedown anchor");
                    }
                }
            }
        }
        if (!probe) std::cout << "PASS: custom upper/lower/interpolated bounds preserve shakedown anchor\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
