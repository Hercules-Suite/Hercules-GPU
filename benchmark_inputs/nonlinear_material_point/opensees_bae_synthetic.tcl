# Standalone OpenSees characterization test for the Hercules vonMises_baE
# material-point fixture. This is a synthetic numerical comparison, not a
# calibrated production soil model.
#
# Usage:
#   OpenSees opensees_bae_synthetic.tcl
#   OpenSees opensees_bae_synthetic.tcl 0.5 opensees_bae_beta05.csv
#
# Optional arguments are the OpenSees integration beta and output CSV path.

set beta 1.0
if {$argc >= 1} {
    set beta [lindex $argv 0]
}
if {![string is double -strict $beta] || $beta < 0.0 || $beta > 1.0} {
    puts stderr "beta must be a number between 0 and 1"
    exit 2
}

set outputPath [format "opensees_bae_beta_%g.csv" $beta]
if {$argc >= 2} {
    set outputPath [lindex $argv 1]
}
if {$argc > 2} {
    puts stderr "usage: OpenSees opensees_bae_synthetic.tcl ?beta? ?output.csv?"
    exit 2
}

# Hercules-to-OpenSees mapping:
#   K = lambda + 2*mu/3, G = mu, Su = cohesion,
#   Ho = 0, h = psi0*G, m = m.
set rho 0.0
set bulkModulus 2.13333333333333343e8
set shearModulus 8.00000000000000000e7
set undrainedShearStrength 2.93938769133981405e4
set initialHardening 0.0
set hardeningScale 7.46480000000000000e7
set hardeningExponent 1.24500000000000000
set k0 0.5

# Engineering gamma_xy corresponding to twice the tensorial epsilon_xy in
# vonmises_bae_opensees_synthetic.in. The repeated zero followed by reloading
# exposes the two implementations' different zero-increment history rules.
set gammaHistory {
     0.000000e+00
     2.500000e-04
     5.000000e-04
     7.500000e-04
     1.000000e-03
     7.500000e-04
     5.000000e-04
     2.500000e-04
     0.000000e+00
    -2.500000e-04
    -5.000000e-04
    -7.500000e-04
    -1.000000e-03
    -7.500000e-04
    -5.000000e-04
    -2.500000e-04
     0.000000e+00
     0.000000e+00
     2.500000e-04
     5.000000e-04
}

wipe
model BasicBuilder -ndm 3 -ndf 3
nDMaterial MultiaxialCyclicPlasticity 1 \
    $rho $bulkModulus $shearModulus $undrainedShearStrength \
    $initialHardening $hardeningScale $hardeningExponent $beta $k0

# Unit cube. The bottom face is fixed; the top face has one shared x
# displacement and fixed y/z displacements, producing affine simple shear.
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

fix 1 1 1 1
fix 2 1 1 1
fix 5 1 1 1
fix 6 1 1 1
fix 3 0 1 1
fix 4 0 1 1
fix 7 0 1 1
fix 8 0 1 1
equalDOF 3 4 1
equalDOF 3 7 1
equalDOF 3 8 1

element stdBrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0
timeSeries Linear 1
pattern Plain 1 1 {
    load 3 1.0 0.0 0.0
}

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-12 50 0
algorithm Newton
integrator LoadControl 0.0
analysis Static

# Prime and commit the material's initial zero-strain state before recording
# the prescribed history.
if {[analyze 1] != 0} {
    puts stderr "OpenSees zero-state initialization failed"
    exit 41
}

set output [open $outputPath w]
puts $output "step,epsilon_xy,gamma_xy,sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_yz,sigma_xz"

set currentGamma 0.0
set step 0
foreach gamma $gammaHistory {
    set increment [expr {$gamma - $currentGamma}]
    if {$increment == 0.0} {
        integrator LoadControl 0.0
    } else {
        integrator DisplacementControl 3 1 $increment
    }

    if {[analyze 1] != 0} {
        puts stderr "OpenSees analysis failed at material step $step"
        close $output
        exit 42
    }

    set stress [eleResponse 1 material 1 stress]
    if {[llength $stress] != 6} {
        puts stderr "OpenSees brick did not return six material stresses"
        close $output
        exit 43
    }
    set epsilonXY [expr {0.5 * $gamma}]
    puts $output [format "%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e" \
        $step $epsilonXY $gamma \
        [lindex $stress 0] [lindex $stress 1] [lindex $stress 2] \
        [lindex $stress 3] [lindex $stress 4] [lindex $stress 5]]

    set currentGamma $gamma
    incr step
}

close $output
puts "OpenSees BAE history complete: $step steps"
puts "beta = $beta"
puts "CSV output: $outputPath"
wipe
exit 0
