# Canonical drained monotonic triaxial-compression characterization for the
# OpenSees Manzari--Dafalias material.  OpenSees uses kN, m, s, and kPa.
#
# Usage:
#   OpenSees opensees_manzari_dafalias_drained_triaxial.tcl
#   OpenSees opensees_manzari_dafalias_drained_triaxial.tcl output.csv

set outputPath "opensees_md_drained_triaxial.csv"
if {$argc >= 1} {
    set outputPath [lindex $argv 0]
}
if {$argc > 1} {
    puts stderr "usage: OpenSees opensees_manzari_dafalias_drained_triaxial.tcl ?output.csv?"
    exit 2
}

set confinementKPa 100.0
set consolidationSteps 50
set axialSteps 150
set finalAxialStrain -0.15

wipe
model BasicBuilder -ndm 3 -ndf 3

# G0 nu e_init Mc c lambda_c e0 ksi P_atm[kPa] m h0 ch nb A0 nd
# z_max cz density integrationScheme tangentType jacobianType TolF TolR
nDMaterial ManzariDafalias 1 \
    125.0 0.05 0.80 1.25 0.712 0.019 0.934 0.70 101.3 \
    0.01 7.05 0.968 1.10 0.704 3.50 4.0 600.0 0.0 \
    1 0 1 1.0e-4 1.0e-7

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

# Symmetry planes remove rigid motion.  The equal-DOF constraints keep each
# loaded face planar while allowing all three normal strains to evolve.
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1
fix 4 1 0 1
fix 5 1 1 0
fix 6 0 1 0
fix 7 0 0 0
fix 8 1 0 0

equalDOF 2 3 1
equalDOF 2 6 1
equalDOF 2 7 1
equalDOF 4 3 2
equalDOF 4 8 2
equalDOF 4 7 2
equalDOF 5 6 3
equalDOF 5 8 3
equalDOF 5 7 3

element stdBrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0
updateMaterialStage -material 1 -stage 0

constraints Transformation
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-11 100 0
algorithm Newton

# A unit cube makes each representative face force numerically equal to the
# applied normal stress.  Compression is negative in OpenSees.
pattern Plain 1 Linear {
    load 2 [expr {-$confinementKPa}] 0.0 0.0
    load 4 0.0 [expr {-$confinementKPa}] 0.0
    load 5 0.0 0.0 [expr {-$confinementKPa}]
}
integrator LoadControl [expr {1.0/$consolidationSteps}]
analysis Static

set output [open $outputPath w]
puts $output "step,phase,eps_xx,eps_yy,eps_zz,eps_xy,eps_yz,eps_xz,sigma_xx_kpa,sigma_yy_kpa,sigma_zz_kpa,sigma_xy_kpa,sigma_yz_kpa,sigma_xz_kpa,p_eff_kpa,q_kpa,q_over_p,alpha_xx,alpha_yy,alpha_zz,alpha_xy,alpha_yz,alpha_xz,fabric_xx,fabric_yy,fabric_zz,fabric_xy,fabric_yz,fabric_xz,alpha_in_xx,alpha_in_yy,alpha_in_zz,alpha_in_xy,alpha_in_yz,alpha_in_xz,void_ratio,plastic_multiplier,confinement_kpa,test_id,driver_version,domain_time"

set step 0
proc writeState {phaseName} {
    global output step confinementKPa
    set stress [eleResponse 1 material 1 stress]
    set state [eleResponse 1 material 1 state]
    if {[llength $stress] != 6 || [llength $state] != 26} {
        puts stderr "OpenSees material response has an unexpected size at step $step"
        close $output
        exit 43
    }
    set exx [nodeDisp 2 1]
    set eyy [nodeDisp 4 2]
    set ezz [nodeDisp 5 3]
    set sxx [lindex $stress 0]
    set syy [lindex $stress 1]
    set szz [lindex $stress 2]
    set sxy [lindex $stress 3]
    set syz [lindex $stress 4]
    set sxz [lindex $stress 5]
    set p [expr {-($sxx+$syy+$szz)/3.0}]
    set dsxx [expr {$sxx+$p}]
    set dsyy [expr {$syy+$p}]
    set dszz [expr {$szz+$p}]
    set q [expr {sqrt(1.5*($dsxx*$dsxx+$dsyy*$dsyy+$dszz*$dszz+2.0*($sxy*$sxy+$syz*$syz+$sxz*$sxz)))}]
    set eta 0.0
    if {$p > 0.0} {
        set eta [expr {$q/$p}]
    }
    puts $output [format "%d,%s,%.17e,%.17e,%.17e,0.0,0.0,0.0,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,md_drained_triaxial,1,%.17e" \
        $step $phaseName $exx $eyy $ezz \
        $sxx $syy $szz $sxy $syz $sxz $p $q $eta \
        [lindex $state 6] [lindex $state 7] [lindex $state 8] \
        [lindex $state 9] [lindex $state 10] [lindex $state 11] \
        [lindex $state 12] [lindex $state 13] [lindex $state 14] \
        [lindex $state 15] [lindex $state 16] [lindex $state 17] \
        [lindex $state 18] [lindex $state 19] [lindex $state 20] \
        [lindex $state 21] [lindex $state 22] [lindex $state 23] \
        [lindex $state 24] [lindex $state 25] $confinementKPa [getTime]]
    incr step
}

writeState initial
for {set i 0} {$i < $consolidationSteps} {incr i} {
    if {[analyze 1] != 0} {
        puts stderr "OpenSees isotropic consolidation failed at increment $i"
        close $output
        exit 42
    }
    writeState consolidation
}

# Freeze the isotropic load.  A dummy axial reference load lets
# DisplacementControl prescribe compression while the x/y face forces remain
# fixed at the confinement pressure.
loadConst -time 0.0
updateMaterialStage -material 1 -stage 1
pattern Plain 2 Linear {
    load 5 0.0 0.0 -1.0
}
integrator DisplacementControl 5 3 [expr {$finalAxialStrain/$axialSteps}]
for {set i 0} {$i < $axialSteps} {incr i} {
    if {[analyze 1] != 0} {
        puts stderr "OpenSees drained triaxial loading failed at increment $i"
        close $output
        exit 44
    }
    writeState monotonic
}

close $output
puts "OpenSees drained Manzari-Dafalias triaxial history complete: $step states"
puts "CSV output: $outputPath"
wipe
exit 0
