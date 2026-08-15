# Standalone OpenSees reference for
# manzari_dafalias_modified_euler_synthetic.in.
#
# Usage:
#   OpenSees opensees_manzari_dafalias_synthetic.tcl
#   OpenSees opensees_manzari_dafalias_synthetic.tcl output.csv
#   OpenSees opensees_manzari_dafalias_synthetic.tcl output.csv forward_euler
#
# OpenSees is run in kPa so the UW low-pressure dilatancy regularization keeps
# its calibrated pressure scale. The companion comparator converts stress to
# Hercules SI units (Pa).

set outputPath "opensees_manzari_dafalias_synthetic.csv"
if {$argc >= 1} {
    set outputPath [lindex $argv 0]
}
set integrationSchemeName "modified_euler"
if {$argc >= 2} {
    set integrationSchemeName [string tolower [lindex $argv 1]]
}
if {$argc > 2} {
    puts stderr "usage: OpenSees opensees_manzari_dafalias_synthetic.tcl ?output.csv? ?modified_euler|forward_euler?"
    exit 2
}
if {$integrationSchemeName eq "modified_euler"} {
    set integrationScheme 1
} elseif {$integrationSchemeName eq "forward_euler"} {
    set integrationScheme 5
} else {
    puts stderr "unknown Manzari-Dafalias integration scheme: $integrationSchemeName"
    exit 2
}

set isotropicStrainHistory {
     0.0000e+00 -2.0000e-04 -4.0000e-04 -6.0000e-04 -8.0000e-04
    -8.0000e-04 -8.0000e-04 -8.0000e-04 -8.0000e-04 -8.0000e-04
    -8.0000e-04 -8.0000e-04 -8.0000e-04 -8.0000e-04 -8.0000e-04
    -8.0000e-04 -8.0000e-04 -8.0000e-04 -8.0000e-04 -8.0000e-04
    -8.0000e-04
}
set gammaXYHistory {
     0.0e+00  0.0e+00  0.0e+00  0.0e+00  0.0e+00
     2.0e-05  4.0e-05  8.0e-05  1.2e-04  1.6e-04
     1.2e-04  8.0e-05  4.0e-05  0.0e+00 -4.0e-05
    -8.0e-05 -1.2e-04 -8.0e-05 -4.0e-05  0.0e+00
     0.0e+00
}
set timeHistory {
     0  1  2  3  4  5  6  7  8  9 10
    11 12 13 14 15 16 17 18 19 20 21
}
set eValues [linsert $isotropicStrainHistory 0 0.0]
set gammaValues [linsert $gammaXYHistory 0 0.0]
set ePlusGammaValues {0.0}
foreach e $isotropicStrainHistory gamma $gammaXYHistory {
    lappend ePlusGammaValues [expr {$e + $gamma}]
}

wipe
model BasicBuilder -ndm 3 -ndf 3

# G0 nu e_init Mc c lambda_c e0 ksi P_atm[kPa] m h0 ch nb A0 nd
# z_max cz density integrationScheme tangentType jacobianType TolF TolR
nDMaterial ManzariDafalias 1 \
    125.0 0.05 0.80 1.25 0.712 0.019 0.934 0.70 101.3 \
    0.01 7.05 0.968 1.10 0.704 3.50 4.0 600.0 0.0 \
    $integrationScheme 0 1 1.0e-4 1.0e-7

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

# Prescribe u_x=e*x+gamma*y, u_y=e*y, and u_z=e*z. This produces the
# fixture's uniform diagonal tensor strain and engineering gamma_xy=2*eps_xy.
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1
fix 4 0 0 1
fix 5 1 1 0
fix 6 0 1 0
fix 7 0 0 0
fix 8 0 0 0

timeSeries Path 1 -time $timeHistory -values $eValues -useLast
pattern Plain 1 1 {
    sp 2 1 1.0
    sp 3 2 1.0
    sp 4 2 1.0
    sp 5 3 1.0
    sp 6 1 1.0
    sp 6 3 1.0
    sp 7 2 1.0
    sp 7 3 1.0
    sp 8 2 1.0
    sp 8 3 1.0
}
timeSeries Path 2 -time $timeHistory -values $gammaValues -useLast
pattern Plain 2 2 {
    sp 4 1 1.0
    sp 8 1 1.0
}
timeSeries Path 3 -time $timeHistory -values $ePlusGammaValues -useLast
pattern Plain 3 3 {
    sp 3 1 1.0
    sp 7 1 1.0
}

element stdBrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0
updateMaterialStage -material 1 -stage 1

constraints Penalty 1.0e18 1.0e18
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-12 50 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set output [open $outputPath w]
puts $output "step,eps_xx,eps_yy,eps_zz,eps_xy,sigma_xx_kpa,sigma_yy_kpa,sigma_zz_kpa,sigma_xy_kpa,sigma_yz_kpa,sigma_xz_kpa,alpha_xx,alpha_yy,alpha_zz,alpha_xy,alpha_yz,alpha_xz,fabric_xx,fabric_yy,fabric_zz,fabric_xy,fabric_yz,fabric_xz,alpha_in_xx,alpha_in_yy,alpha_in_zz,alpha_in_xy,alpha_in_yz,alpha_in_xz,void_ratio,plastic_multiplier,integration_scheme_id,driver_version,domain_time"

set step 0
foreach e $isotropicStrainHistory gamma $gammaXYHistory {
    if {[analyze 1] != 0} {
        puts stderr "OpenSees analysis failed at material step $step"
        close $output
        exit 42
    }
    set expectedEPlusGamma [expr {$e + $gamma}]
    if {abs([nodeDisp 2 1] - $e) > 1.0e-12 ||
        abs([nodeDisp 4 1] - $gamma) > 1.0e-12 ||
        abs([nodeDisp 3 1] - $expectedEPlusGamma) > 1.0e-12 ||
        abs([nodeDisp 3 2] - $e) > 1.0e-12 ||
        abs([nodeDisp 5 3] - $e) > 1.0e-12} {
        puts stderr "OpenSees affine displacement mismatch at material step $step"
        close $output
        exit 44
    }
    set stress [eleResponse 1 material 1 stress]
    set state [eleResponse 1 material 1 state]
    if {[llength $stress] != 6 || [llength $state] != 26} {
        puts stderr "OpenSees material response has an unexpected size at step $step"
        close $output
        exit 43
    }
    puts $output [format "%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%d,%.17e" \
        $step $e $e $e [expr {0.5*$gamma}] \
        [lindex $stress 0] [lindex $stress 1] [lindex $stress 2] \
        [lindex $stress 3] [lindex $stress 4] [lindex $stress 5] \
        [lindex $state 6] [lindex $state 7] [lindex $state 8] \
        [lindex $state 9] [lindex $state 10] [lindex $state 11] \
        [lindex $state 12] [lindex $state 13] [lindex $state 14] \
        [lindex $state 15] [lindex $state 16] [lindex $state 17] \
        [lindex $state 18] [lindex $state 19] [lindex $state 20] \
        [lindex $state 21] [lindex $state 22] [lindex $state 23] \
        [lindex $state 24] [lindex $state 25] \
        $integrationScheme 2 [getTime]]
    incr step
}

close $output
puts "OpenSees Manzari-Dafalias $integrationSchemeName history complete: $step steps"
puts "CSV output: $outputPath"
wipe
exit 0
