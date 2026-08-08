# Coupled, sealed, one-element undrained triaxial characterization using the
# OpenSees SSPbrickUP element and Manzari--Dafalias effective-stress material.
# OpenSees units are kN, m, s, kPa, and tonne/m^3.
#
# Usage:
#   OpenSees opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl
#   OpenSees opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl output.csv

set outputPath "opensees_md_sspbrickup_undrained_triaxial.csv"
if {$argc >= 1} {
    set outputPath [lindex $argv 0]
}
if {$argc > 1} {
    puts stderr "usage: OpenSees opensees_manzari_dafalias_sspbrickup_undrained_triaxial.tcl ?output.csv?"
    exit 2
}

set confinementKPa 100.0
set consolidationSteps 50
set axialSteps 200
set finalAxialIncrement -0.05
set maxSubstepDepth 10
set fluidBulkKPa 2.2e6
set fluidDensity 1.0
set voidRatio 0.80
set porosity [expr {$voidRatio/(1.0+$voidRatio)}]
# SSPbrickUP assumes unit Biot coupling. Its input named alpha is instead the
# pressure-stabilization coefficient used by the single-point element.
set biotAlpha 1.0
set pressureStabilization 1.5e-9
set storageModulusKPa [expr {$fluidBulkKPa/$porosity}]
set pressureUniformToleranceKPa 1.0e-6
set storageResidualToleranceKPa 1.0e-3
set confinementToleranceKPa 1.0e-2
set initialUndrainedVolumetricStrain 0.0
set initialUndrainedPressureKPa 0.0

wipe
model BasicBuilder -ndm 3 -ndf 4

# G0 nu e_init Mc c lambda_c e0 ksi P_atm[kPa] m h0 ch nb A0 nd
# z_max cz density integrationScheme tangentType jacobianType TolF TolR
nDMaterial ManzariDafalias 1 \
    125.0 0.05 $voidRatio 1.25 0.712 0.019 0.934 0.70 101.3 \
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

# Symmetry constraints for a conventional triaxial specimen. During drained
# consolidation all eight pore-pressure DOFs are fixed to zero; removing those
# SP constraints seals the element for the coupled undrained stage. Distributed
# face tractions and identical top-face axial displacements produce the affine
# one-element mode without mixing an SP constraint with an equalDOF retained
# node (a combination unsupported by OpenSees' Transformation handler).
fix 1 1 1 1 1
fix 2 0 1 1 1
fix 3 0 0 1 1
fix 4 1 0 1 1
fix 5 1 1 0 1
fix 6 0 1 0 1
fix 7 0 0 0 1
fix 8 1 0 0 1

# Zero permeability gives a sealed element after the pressure SP constraints
# are removed. -lumped row-sums the SSP pressure capacity to match Hercules'
# nodal BrickUP storage volume. Pressure stabilization has zero row sum for
# this uniform mode.
element SSPbrickUP 1 1 2 3 4 5 6 7 8 1 \
    $fluidBulkKPa $fluidDensity 0.0 0.0 0.0 \
    $voidRatio $pressureStabilization -lumped
updateMaterialStage -material 1 -stage 0

constraints Transformation
numberer RCM
system BandGeneral
test NormUnbalance 1.0e-8 100 0
algorithm Newton
integrator LoadControl [expr {1.0/$consolidationSteps}]
analysis Static

set pNode [expr {-$confinementKPa/4.0}]
pattern Plain 1 Linear {
    load 2 $pNode 0.0 0.0 0.0
    load 3 $pNode $pNode 0.0 0.0
    load 4 0.0 $pNode 0.0 0.0
    load 5 0.0 0.0 $pNode 0.0
    load 6 $pNode 0.0 $pNode 0.0
    load 7 $pNode $pNode $pNode 0.0
    load 8 0.0 $pNode $pNode 0.0
}

set output [open $outputPath w]
puts $output "step,phase,eps_xx,eps_yy,eps_zz,eps_xy,eps_yz,eps_xz,sigma_xx_kpa,sigma_yy_kpa,sigma_zz_kpa,sigma_xy_kpa,sigma_yz_kpa,sigma_xz_kpa,pore_pressure_kpa,total_sigma_xx_kpa,total_sigma_yy_kpa,total_sigma_zz_kpa,total_sigma_xy_kpa,total_sigma_yz_kpa,total_sigma_xz_kpa,p_eff_kpa,q_kpa,q_over_p,ru,alpha_xx,alpha_yy,alpha_zz,alpha_xy,alpha_yz,alpha_xz,fabric_xx,fabric_yy,fabric_zz,fabric_xy,fabric_yz,fabric_xz,alpha_in_xx,alpha_in_yy,alpha_in_zz,alpha_in_xy,alpha_in_yz,alpha_in_xz,void_ratio,plastic_multiplier,confinement_kpa,fluid_bulk_modulus_kpa,porosity,biot_alpha,pressure_stabilization,test_id,driver_version,domain_time"

set step 0
proc writeState {phaseName} {
    global output step confinementKPa fluidBulkKPa porosity biotAlpha
    global pressureStabilization
    global storageModulusKPa pressureUniformToleranceKPa
    global storageResidualToleranceKPa confinementToleranceKPa
    global initialUndrainedVolumetricStrain initialUndrainedPressureKPa

    set stress [eleResponse 1 stress3D6]
    set strain [eleResponse 1 strain3D6]
    set state [eleResponse 1 state]
    set alpha [eleResponse 1 alpha]
    set fabric [eleResponse 1 fabric]
    set alphaIn [eleResponse 1 alpha_in]
    if {[llength $stress] != 6 || [llength $strain] != 6 ||
        [llength $state] != 26 || [llength $alpha] != 6 ||
        [llength $fabric] != 6 || [llength $alphaIn] != 6} {
        puts stderr "OpenSees SSPbrickUP material response has an unexpected size at step $step"
        close $output
        exit 51
    }

    set exx [lindex $strain 0]
    set eyy [lindex $strain 1]
    set ezz [lindex $strain 2]
    set exy [expr {0.5*[lindex $strain 3]}]
    set eyz [expr {0.5*[lindex $strain 4]}]
    set exz [expr {0.5*[lindex $strain 5]}]
    set affineError [expr {max(abs($exy), abs($eyz), abs($exz))}]
    foreach node {2 3 6 7} {
        set affineError [expr {max($affineError,
                                   abs([nodeDisp $node 1]-$exx))}]
    }
    foreach node {3 4 7 8} {
        set affineError [expr {max($affineError,
                                   abs([nodeDisp $node 2]-$eyy))}]
    }
    foreach node {5 6 7 8} {
        set affineError [expr {max($affineError,
                                   abs([nodeDisp $node 3]-$ezz))}]
    }
    if {$affineError > 1.0e-12} {
        puts stderr "OpenSees SSPbrickUP deformation is not affine at step $step: error=$affineError"
        close $output
        exit 52
    }

    set pressureValues {}
    for {set node 1} {$node <= 8} {incr node} {
        lappend pressureValues [nodeVel $node 4]
    }
    set pressure 0.0
    foreach value $pressureValues {
        set pressure [expr {$pressure+$value/8.0}]
    }
    set pressureSpread 0.0
    foreach value $pressureValues {
        set pressureSpread [expr {max($pressureSpread,abs($value-$pressure))}]
    }
    if {$pressureSpread > $pressureUniformToleranceKPa} {
        puts stderr "OpenSees SSPbrickUP pressure is not uniform at step $step: spread=$pressureSpread kPa"
        close $output
        exit 53
    }

    set sxx [lindex $stress 0]
    set syy [lindex $stress 1]
    set szz [lindex $stress 2]
    set sxy [lindex $stress 3]
    set syz [lindex $stress 4]
    set sxz [lindex $stress 5]
    set txx [expr {$sxx-$biotAlpha*$pressure}]
    set tyy [expr {$syy-$biotAlpha*$pressure}]
    set tzz [expr {$szz-$biotAlpha*$pressure}]
    set txy $sxy
    set tyz $syz
    set txz $sxz

    if {$phaseName eq "coupled"} {
        set volumetricStrain [expr {$exx+$eyy+$ezz}]
        set expectedPressure [expr {$initialUndrainedPressureKPa-
            $storageModulusKPa*($volumetricStrain-
                                $initialUndrainedVolumetricStrain)}]
        set storageResidual [expr {abs($pressure-$expectedPressure)}]
        if {$storageResidual > $storageResidualToleranceKPa} {
            puts stderr "OpenSees SSPbrickUP sealed-storage residual exceeds $storageResidualToleranceKPa kPa at step $step: error=$storageResidual"
            close $output
            exit 54
        }
        set confinementError [expr {max(abs(-$txx-$confinementKPa),
                                        abs(-$tyy-$confinementKPa))}]
        if {$confinementError > $confinementToleranceKPa} {
            puts stderr "OpenSees SSPbrickUP total radial confinement residual exceeds $confinementToleranceKPa kPa at step $step: error=$confinementError"
            close $output
            exit 55
        }
    }

    set pEff [expr {-($sxx+$syy+$szz)/3.0}]
    set dsxx [expr {$sxx+$pEff}]
    set dsyy [expr {$syy+$pEff}]
    set dszz [expr {$szz+$pEff}]
    set q [expr {sqrt(1.5*($dsxx*$dsxx+$dsyy*$dsyy+$dszz*$dszz+
                           2.0*($sxy*$sxy+$syz*$syz+$sxz*$sxz)))}]
    set eta 0.0
    if {$pEff > 0.0} {
        set eta [expr {$q/$pEff}]
    }
    set ru [expr {$pressure/$confinementKPa}]

    set row [list $step $phaseName]
    set numeric [concat \
        [list $exx $eyy $ezz $exy $eyz $exz \
              $sxx $syy $szz $sxy $syz $sxz $pressure \
              $txx $tyy $tzz $txy $tyz $txz \
              $pEff $q $eta $ru] \
        $alpha $fabric $alphaIn \
        [list [lindex $state 24] [lindex $state 25] \
              $confinementKPa $fluidBulkKPa $porosity $biotAlpha \
              $pressureStabilization]]
    foreach value $numeric {
        lappend row [format "%.17e" $value]
    }
    lappend row md_sspbrickup_undrained_triaxial 2 [format "%.17e" [getTime]]
    puts $output [join $row ,]
    incr step
}

writeState initial
for {set i 0} {$i < $consolidationSteps} {incr i} {
    if {[analyze 1] != 0} {
        puts stderr "OpenSees SSPbrickUP consolidation failed at increment $i"
        close $output
        exit 56
    }
    writeState consolidation
}

# Freeze the total confinement loads, open the pressure DOFs from their
# drained-zero consolidation constraints, and enable plasticity.
loadConst -time 0.0
for {set node 1} {$node <= 8} {incr node} {
    remove sp $node 4
}
updateMaterialStage -material 1 -stage 1

set boundaryStrain [eleResponse 1 strain3D6]
set initialUndrainedVolumetricStrain [expr {
    [lindex $boundaryStrain 0]+[lindex $boundaryStrain 1]+
    [lindex $boundaryStrain 2]}]
set initialUndrainedPressureKPa 0.0
set initialAxialDisplacement [nodeDisp 7 3]
set finalAxialDisplacement [expr {$initialAxialDisplacement+
                                  $finalAxialIncrement}]

# The Path series prescribes the absolute axial displacement, beginning at the
# consolidated value so stage activation introduces no displacement jump.
timeSeries Path 2 -time {0.0 1.0} \
    -values [list $initialAxialDisplacement $finalAxialDisplacement] -useLast
pattern Plain 2 2 {
    sp 5 3 1.0
    sp 6 3 1.0
    sp 7 3 1.0
    sp 8 3 1.0
}

# Rebuild the analysis after changing the pressure and axial SP constraints.
# The model, committed material state, constant confinement loads, and newly
# defined displacement pattern remain in the Domain across wipeAnalysis.
wipeAnalysis
constraints Transformation
numberer RCM
system BandGeneral
test NormUnbalance 1.0e-3 100 0
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set acceptedCoupledSubsteps 0
proc advanceCoupledIncrement {increment depth} {
    global maxSubstepDepth acceptedCoupledSubsteps
    test NormUnbalance 1.0e-3 100 0
    algorithm Newton
    set ok [analyze 1 $increment]
    if {$ok == 0} {
        writeState coupled
        incr acceptedCoupledSubsteps
        return 0
    }
    if {$depth >= $maxSubstepDepth} {
        puts stderr "OpenSees SSPbrickUP increment failed at maximum substep depth $depth"
        return -10
    }
    set halfIncrement [expr {0.5*$increment}]
    puts stderr "OpenSees SSPbrickUP increment did not converge; subdividing to $halfIncrement (depth [expr {$depth+1}])"
    set ok [advanceCoupledIncrement $halfIncrement [expr {$depth+1}]]
    if {$ok != 0} {
        return $ok
    }
    return [advanceCoupledIncrement $halfIncrement [expr {$depth+1}]]
}

set nominalDt [expr {1.0/$axialSteps}]
for {set i 0} {$i < $axialSteps} {incr i} {
    set ok [advanceCoupledIncrement $nominalDt 0]
    if {$ok != 0} {
        puts stderr "OpenSees SSPbrickUP undrained loading failed at increment $i"
        close $output
        exit 57
    }
}
writeState complete

close $output
puts "OpenSees SSPbrickUP Manzari-Dafalias history complete: $step states"
puts "Accepted coupled increments (including substeps): $acceptedCoupledSubsteps"
puts "CSV output: $outputPath"
wipe
exit 0
