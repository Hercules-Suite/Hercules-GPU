# Canonical constant-volume cyclic-simple-shear characterization for the
# OpenSees Manzari--Dafalias material.  OpenSees uses kN, m, s, and kPa.
# Excess pore pressure is the pressure required to maintain constant total
# mean stress: delta_u = p'_0 - p'.  This is a material-level undrained test,
# not a u-p finite-element benchmark.
#
# Usage:
#   OpenSees opensees_manzari_dafalias_undrained_cyclic_shear.tcl
#   OpenSees opensees_manzari_dafalias_undrained_cyclic_shear.tcl output.csv

set outputPath "opensees_md_undrained_cyclic_shear.csv"
if {$argc >= 1} {
    set outputPath [lindex $argv 0]
}
if {$argc > 1} {
    puts stderr "usage: OpenSees opensees_manzari_dafalias_undrained_cyclic_shear.tcl ?output.csv?"
    exit 2
}

set consolidationSteps 26
set consolidationStrain -2.6e-3
set cyclicAmplitude 5.0e-3
set cyclicIncrement 2.5e-4
set numberOfCycles 5

# Build the total-strain path.  The diagonal strain is held fixed during the
# shear phase, so every cyclic increment has zero volumetric strain.
set normalHistory {0.0}
set gammaHistory {0.0}
set phaseHistory {initial}
for {set i 1} {$i <= $consolidationSteps} {incr i} {
    lappend normalHistory [expr {$consolidationStrain*$i/$consolidationSteps}]
    lappend gammaHistory 0.0
    lappend phaseHistory consolidation
}

# The repeated zero-shear state is the first stage-1 update.
lappend normalHistory $consolidationStrain
lappend gammaHistory 0.0
lappend phaseHistory cyclic

set turningPoints {}
for {set cycle 0} {$cycle < $numberOfCycles} {incr cycle} {
    lappend turningPoints $cyclicAmplitude [expr {-$cyclicAmplitude}]
}
lappend turningPoints 0.0
set currentGamma 0.0
foreach targetGamma $turningPoints {
    set segmentSteps [expr {int(round(abs($targetGamma-$currentGamma)/$cyclicIncrement))}]
    for {set i 1} {$i <= $segmentSteps} {incr i} {
        lappend normalHistory $consolidationStrain
        lappend gammaHistory [expr {$currentGamma+($targetGamma-$currentGamma)*$i/$segmentSteps}]
        lappend phaseHistory cyclic
    }
    set currentGamma $targetGamma
}

set numberOfStates [llength $normalHistory]
set timeHistory {}
set ePlusGammaHistory {}
for {set i 0} {$i < $numberOfStates} {incr i} {
    lappend timeHistory $i
    lappend ePlusGammaHistory [expr {[lindex $normalHistory $i]+[lindex $gammaHistory $i]}]
}

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

# Prescribe u_x=e*x+gamma*y, u_y=e*y, and u_z=e*z.
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1
fix 4 0 0 1
fix 5 1 1 0
fix 6 0 1 0
fix 7 0 0 0
fix 8 0 0 0

timeSeries Path 1 -time $timeHistory -values $normalHistory -useLast
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
timeSeries Path 2 -time $timeHistory -values $gammaHistory -useLast
pattern Plain 2 2 {
    sp 4 1 1.0
    sp 8 1 1.0
}
timeSeries Path 3 -time $timeHistory -values $ePlusGammaHistory -useLast
pattern Plain 3 3 {
    sp 3 1 1.0
    sp 7 1 1.0
}

element stdBrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0
updateMaterialStage -material 1 -stage 0

constraints Penalty 1.0e18 1.0e18
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-12 100 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

set output [open $outputPath w]
puts $output "step,phase,eps_xx,eps_yy,eps_zz,eps_xy,eps_yz,eps_xz,sigma_xx_kpa,sigma_yy_kpa,sigma_zz_kpa,sigma_xy_kpa,sigma_yz_kpa,sigma_xz_kpa,p_eff_kpa,q_kpa,excess_pore_pressure_kpa,ru,alpha_xx,alpha_yy,alpha_zz,alpha_xy,alpha_yz,alpha_xz,fabric_xx,fabric_yy,fabric_zz,fabric_xy,fabric_yz,fabric_xz,alpha_in_xx,alpha_in_yy,alpha_in_zz,alpha_in_xy,alpha_in_yz,alpha_in_xz,void_ratio,plastic_multiplier,test_id,driver_version,domain_time"

set step 0
set initialEffectivePressure 0.0
set pressureInitialized 0
proc writeState {phaseName} {
    global output step initialEffectivePressure pressureInitialized
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
    set exy [expr {0.5*[nodeDisp 4 1]}]
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
    if {$phaseName eq "cyclic" && !$pressureInitialized} {
        set initialEffectivePressure $p
        set pressureInitialized 1
    }
    set excessPressure 0.0
    set ru 0.0
    if {$pressureInitialized} {
        set excessPressure [expr {$initialEffectivePressure-$p}]
        if {$initialEffectivePressure > 0.0} {
            set ru [expr {$excessPressure/$initialEffectivePressure}]
        }
    }
    puts $output [format "%d,%s,%.17e,%.17e,%.17e,%.17e,0.0,0.0,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,md_undrained_cyclic_shear,1,%.17e" \
        $step $phaseName $exx $eyy $ezz $exy \
        $sxx $syy $szz $sxy $syz $sxz $p $q $excessPressure $ru \
        [lindex $state 6] [lindex $state 7] [lindex $state 8] \
        [lindex $state 9] [lindex $state 10] [lindex $state 11] \
        [lindex $state 12] [lindex $state 13] [lindex $state 14] \
        [lindex $state 15] [lindex $state 16] [lindex $state 17] \
        [lindex $state 18] [lindex $state 19] [lindex $state 20] \
        [lindex $state 21] [lindex $state 22] [lindex $state 23] \
        [lindex $state 24] [lindex $state 25] [getTime]]
    incr step
}

writeState initial
for {set index 1} {$index < $numberOfStates} {incr index} {
    if {$index == [expr {$consolidationSteps+1}]} {
        updateMaterialStage -material 1 -stage 1
    }
    if {[analyze 1] != 0} {
        puts stderr "OpenSees undrained cyclic analysis failed at state $index"
        close $output
        exit 42
    }
    writeState [lindex $phaseHistory $index]
}

close $output
puts "OpenSees undrained Manzari-Dafalias cyclic history complete: $step states"
puts "CSV output: $outputPath"
wipe
exit 0
