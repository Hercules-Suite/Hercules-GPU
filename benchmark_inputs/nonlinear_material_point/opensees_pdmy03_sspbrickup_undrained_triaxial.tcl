# Coupled, sealed, one-element undrained triaxial characterization using
# OpenSees SSPbrickUP and PressureDependMultiYield03. Units are kN, m, s,
# kPa, and tonne/m^3.
#
# Usage:
#   OpenSees opensees_pdmy03_sspbrickup_undrained_triaxial.tcl ?output.csv?

set outputPath "opensees_pdmy03_sspbrickup_undrained_triaxial.csv"
if {$argc >= 1} { set outputPath [lindex $argv 0] }
if {$argc > 1} {
    puts stderr "usage: OpenSees opensees_pdmy03_sspbrickup_undrained_triaxial.tcl ?output.csv?"
    exit 2
}

set confinementKPa 100.0
set consolidationSteps 50
set axialSteps 240
set axialTimes {0.0 0.15 0.30 0.45 0.60 0.75 1.0}
set axialIncrements {0.0 -0.005 0.005 -0.010 0.010 -0.020 -0.040}
set maxSubstepDepth 10
set fluidBulkKPa 2.2e6
set fluidDensity 1.0
set voidRatio 0.80
set porosity [expr {$voidRatio/(1.0+$voidRatio)}]
set biotAlpha 1.0
set pressureStabilization 1.5e-9
set storageModulusKPa [expr {$fluidBulkKPa/$porosity}]
set pressureUniformToleranceKPa 1.0e-6
set storageResidualToleranceKPa 2.0e-3
set confinementToleranceKPa 2.0e-2
set affineTolerance 1.0e-8
set initialUndrainedVolumetricStrain 0.0
set initialUndrainedPressureKPa 0.0

wipe
model BasicBuilder -ndm 3 -ndf 4

# nd rho Gref Kref phi peakShearStrain refPressure pressureDependence
# phaseTransform contraction1 contraction3 dilation1 dilation3 liquefaction1
# liquefaction2 liquefaction4 numberOfSurfaces cohesion pressureDependence2
# atmosphericPressure density
nDMaterial PressureDependMultiYield03 1 3 0.0 \
    9.0e4 2.2e5 32.0 0.10 101.3 0.50 26.0 0.067 0.50 \
    0.23 1.0 0.0 0.06 3.0 0.27 20 10.0 0.01 101.3 1.73

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

fix 1 1 1 1 1
fix 2 0 1 1 1
fix 3 0 0 1 1
fix 4 1 0 1 1
fix 5 1 1 0 1
fix 6 0 1 0 1
fix 7 0 0 0 1
fix 8 1 0 0 1

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
puts $output "step,phase,eps_xx,eps_yy,eps_zz,eps_xy,eps_yz,eps_xz,sigma_xx_kpa,sigma_yy_kpa,sigma_zz_kpa,sigma_xy_kpa,sigma_yz_kpa,sigma_xz_kpa,pore_pressure_kpa,total_sigma_xx_kpa,total_sigma_yy_kpa,total_sigma_zz_kpa,total_sigma_xy_kpa,total_sigma_yz_kpa,total_sigma_xz_kpa,p_eff_kpa,q_kpa,q_over_p,ru,confinement_kpa,fluid_bulk_modulus_kpa,porosity,biot_alpha,pressure_stabilization,test_id,driver_version,domain_time"

set step 0
proc writeState {phaseName} {
    global output step confinementKPa fluidBulkKPa porosity biotAlpha
    global pressureStabilization storageModulusKPa
    global pressureUniformToleranceKPa storageResidualToleranceKPa
    global confinementToleranceKPa affineTolerance
    global initialUndrainedVolumetricStrain initialUndrainedPressureKPa

    set stress [eleResponse 1 stress3D6]
    set strain [eleResponse 1 strain3D6]
    if {[llength $stress] < 6 || [llength $strain] != 6} {
        puts stderr "OpenSees PDMY03 SSPbrickUP response has an unexpected size at step $step"
        close $output
        exit 51
    }

    set exx [lindex $strain 0]
    set eyy [lindex $strain 1]
    set ezz [lindex $strain 2]
    set exy [expr {0.5*[lindex $strain 3]}]
    set eyz [expr {0.5*[lindex $strain 4]}]
    set exz [expr {0.5*[lindex $strain 5]}]
    set affineError [expr {max(abs($exy),abs($eyz),abs($exz))}]
    foreach node {2 3 6 7} {
        set affineError [expr {max($affineError,abs([nodeDisp $node 1]-$exx))}]
    }
    foreach node {3 4 7 8} {
        set affineError [expr {max($affineError,abs([nodeDisp $node 2]-$eyy))}]
    }
    foreach node {5 6 7 8} {
        set affineError [expr {max($affineError,abs([nodeDisp $node 3]-$ezz))}]
    }
    if {$affineError > $affineTolerance} {
        puts stderr "OpenSees PDMY03 SSPbrickUP deformation is not affine at step $step: error=$affineError"
        close $output
        exit 52
    }

    set pressure 0.0
    set pressureValues {}
    for {set node 1} {$node <= 8} {incr node} {
        set value [nodeVel $node 4]
        lappend pressureValues $value
        set pressure [expr {$pressure+$value/8.0}]
    }
    set pressureSpread 0.0
    foreach value $pressureValues {
        set pressureSpread [expr {max($pressureSpread,abs($value-$pressure))}]
    }
    if {$pressureSpread > $pressureUniformToleranceKPa} {
        puts stderr "OpenSees PDMY03 SSPbrickUP pressure is not uniform at step $step: spread=$pressureSpread kPa"
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

    if {$phaseName eq "coupled"} {
        set volumetricStrain [expr {$exx+$eyy+$ezz}]
        set expectedPressure [expr {$initialUndrainedPressureKPa-
            $storageModulusKPa*($volumetricStrain-
                                $initialUndrainedVolumetricStrain)}]
        set storageResidual [expr {abs($pressure-$expectedPressure)}]
        if {$storageResidual > $storageResidualToleranceKPa} {
            puts stderr "OpenSees PDMY03 sealed-storage residual exceeds $storageResidualToleranceKPa kPa at step $step: error=$storageResidual"
            close $output
            exit 54
        }
        set confinementError [expr {max(abs(-$txx-$confinementKPa),
                                        abs(-$tyy-$confinementKPa))}]
        if {$confinementError > $confinementToleranceKPa} {
            puts stderr "OpenSees PDMY03 total radial confinement residual exceeds $confinementToleranceKPa kPa at step $step: error=$confinementError"
            close $output
            exit 55
        }
    }

    set pEff [expr {-($sxx+$syy+$szz)/3.0}]
    set dxx [expr {$sxx+$pEff}]
    set dyy [expr {$syy+$pEff}]
    set dzz [expr {$szz+$pEff}]
    set q [expr {sqrt(1.5*($dxx*$dxx+$dyy*$dyy+$dzz*$dzz+
                           2.0*($sxy*$sxy+$syz*$syz+$sxz*$sxz)))}]
    set eta [expr {$pEff > 0.0 ? $q/$pEff : 0.0}]
    set ru [expr {$pressure/$confinementKPa}]
    set row [list $step $phaseName]
    foreach value [list $exx $eyy $ezz $exy $eyz $exz \
                        $sxx $syy $szz $sxy $syz $sxz $pressure \
                        $txx $tyy $tzz $sxy $syz $sxz \
                        $pEff $q $eta $ru $confinementKPa $fluidBulkKPa \
                        $porosity $biotAlpha $pressureStabilization] {
        lappend row [format "%.17e" $value]
    }
    lappend row pdmy03_sspbrickup_undrained_triaxial 1 \
        [format "%.17e" [getTime]]
    puts $output [join $row ,]
    incr step
}

writeState initial
for {set i 0} {$i < $consolidationSteps} {incr i} {
    if {[analyze 1] != 0} {
        puts stderr "OpenSees PDMY03 SSPbrickUP consolidation failed at increment $i"
        close $output
        exit 56
    }
    writeState consolidation
}

loadConst -time 0.0
for {set node 1} {$node <= 8} {incr node} { remove sp $node 4 }
updateMaterialStage -material 1 -stage 1
set boundaryStrain [eleResponse 1 strain3D6]
set initialUndrainedVolumetricStrain [expr {
    [lindex $boundaryStrain 0]+[lindex $boundaryStrain 1]+
    [lindex $boundaryStrain 2]}]
set initialUndrainedPressureKPa 0.0
set initialAxialDisplacement [nodeDisp 7 3]
set axialValues {}
foreach increment $axialIncrements {
    lappend axialValues [expr {$initialAxialDisplacement+$increment}]
}
timeSeries Path 2 -time $axialTimes -values $axialValues -useLast
pattern Plain 2 2 {
    sp 5 3 1.0
    sp 6 3 1.0
    sp 7 3 1.0
    sp 8 3 1.0
}

wipeAnalysis
constraints Penalty 1.0e18 1.0e18
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-9 100 0
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set acceptedCoupledSubsteps 0
proc advanceCoupledIncrement {increment depth} {
    global maxSubstepDepth acceptedCoupledSubsteps
    test NormDispIncr 1.0e-9 100 0
    algorithm Newton
    set ok [analyze 1 $increment]
    if {$ok == 0} {
        writeState coupled
        incr acceptedCoupledSubsteps
        return 0
    }
    if {$depth >= $maxSubstepDepth} { return -10 }
    set halfIncrement [expr {0.5*$increment}]
    set ok [advanceCoupledIncrement $halfIncrement [expr {$depth+1}]]
    if {$ok != 0} { return $ok }
    return [advanceCoupledIncrement $halfIncrement [expr {$depth+1}]]
}

set nominalDt [expr {1.0/$axialSteps}]
for {set i 0} {$i < $axialSteps} {incr i} {
    if {[advanceCoupledIncrement $nominalDt 0] != 0} {
        puts stderr "OpenSees PDMY03 SSPbrickUP loading failed at increment $i"
        close $output
        exit 57
    }
}
writeState complete
close $output
puts "OpenSees SSPbrickUP PDMY03 undrained triaxial complete: $step states"
puts "Accepted coupled increments (including substeps): $acceptedCoupledSubsteps"
puts "CSV output: $outputPath"
wipe
exit 0
