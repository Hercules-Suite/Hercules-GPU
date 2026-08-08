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
set axialSteps 300
set finalAxialStrain -0.15
set maxSubstepDepth 10

wipe
model BasicBuilder -ndm 3 -ndf 3

# G0 nu e_init Mc c lambda_c e0 ksi P_atm[kPa] m h0 ch nb A0 nd
# z_max cz density integrationScheme tangentType jacobianType TolF TolR
# Integration scheme 1 is Modified Euler.  Tangent type 2 returns the
# accumulated consistent elastoplastic tangent; the elastic tangent (type 0)
# is inadequate for the lateral stress-control solve after yielding begins.
nDMaterial ManzariDafalias 1 \
    125.0 0.05 0.80 1.25 0.712 0.019 0.934 0.70 101.3 \
    0.01 7.05 0.968 1.10 0.704 3.50 4.0 600.0 0.0 \
    1 2 1 1.0e-4 1.0e-7

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

# Symmetry planes remove rigid motion.  Node 7 is the one retained control
# node for all three normal face displacements.  Each other corner appears in
# exactly one equalDOF command; using several retained nodes for different
# DOFs of the same constrained corner does not define a valid affine MP map.
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1
fix 4 1 0 1
fix 5 1 1 0
fix 6 0 1 0
fix 7 0 0 0
fix 8 1 0 0

equalDOF 7 2 1
equalDOF 7 3 1 2
equalDOF 7 4 2
equalDOF 7 5 3
equalDOF 7 6 1 3
equalDOF 7 8 2 3

# The imposed affine normal deformation is uniform, so stabilized single-point
# integration evaluates the same material path as the full-integration brick
# without repeating it at eight identical Gauss points.
element SSPbrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0
updateMaterialStage -material 1 -stage 0

constraints Transformation
numberer RCM
system BandGeneral
# The elastic consolidation stage readily supports a tight global tolerance.
test NormDispIncr 1.0e-10 100 0
algorithm Newton

# A unit cube makes each representative face force numerically equal to the
# applied normal stress.  Compression is negative in OpenSees.
pattern Plain 1 Linear {
    load 7 [expr {-$confinementKPa}] \
           [expr {-$confinementKPa}] \
           [expr {-$confinementKPa}]
}
integrator LoadControl [expr {1.0/$consolidationSteps}]
analysis Static

set output [open $outputPath w]
puts $output "step,phase,eps_xx,eps_yy,eps_zz,eps_xy,eps_yz,eps_xz,sigma_xx_kpa,sigma_yy_kpa,sigma_zz_kpa,sigma_xy_kpa,sigma_yz_kpa,sigma_xz_kpa,p_eff_kpa,q_kpa,q_over_p,alpha_xx,alpha_yy,alpha_zz,alpha_xy,alpha_yz,alpha_xz,fabric_xx,fabric_yy,fabric_zz,fabric_xy,fabric_yz,fabric_xz,alpha_in_xx,alpha_in_yy,alpha_in_zz,alpha_in_xy,alpha_in_yz,alpha_in_xz,void_ratio,plastic_multiplier,confinement_kpa,test_id,driver_version,domain_time"

set step 0
proc writeState {phaseName} {
    global output step confinementKPa
    # SSPbrick delegates direct recorder arguments to its one nDMaterial.  The
    # explicit 3D6 names return the actual center-point tensor responses; unlike
    # stdBrick, SSPbrick does not accept the "material 1" prefix.
    set stress [eleResponse 1 stress3D6]
    set strain [eleResponse 1 strain3D6]
    set state [eleResponse 1 state]
    set alpha [eleResponse 1 alpha]
    set fabric [eleResponse 1 fabric]
    set alphaIn [eleResponse 1 alpha_in]
    if {[llength $stress] != 6 || [llength $strain] != 6 ||
        [llength $state] != 26 || [llength $alpha] != 6 ||
        [llength $fabric] != 6 || [llength $alphaIn] != 6} {
        puts stderr "OpenSees material response has an unexpected size at step $step"
        close $output
        exit 43
    }
    set exx [lindex $strain 0]
    set eyy [lindex $strain 1]
    set ezz [lindex $strain 2]
    # OpenSees stores engineering shear strain; Hercules stores tensor shear.
    set exy [expr {0.5*[lindex $strain 3]}]
    set eyz [expr {0.5*[lindex $strain 4]}]
    set exz [expr {0.5*[lindex $strain 5]}]
    set affineError [expr {max(abs($exx-[nodeDisp 7 1]),
                               abs($eyy-[nodeDisp 7 2]),
                               abs($ezz-[nodeDisp 7 3]),
                               abs($exy), abs($eyz), abs($exz))}]
    if {$affineError > 1.0e-12} {
        puts stderr "OpenSees drained deformation is not affine at step $step: error=$affineError"
        close $output
        exit 45
    }
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
    puts $output [format "%d,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,md_drained_triaxial,7,%.17e" \
        $step $phaseName $exx $eyy $ezz $exy $eyz $exz \
        $sxx $syy $szz $sxy $syz $sxz $p $q $eta \
        {*}$alpha {*}$fabric {*}$alphaIn \
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
# Match the 1.0e-5 global displacement-increment tolerance used by the
# published OpenSees Manzari--Dafalias triaxial example.  This is an outer
# equilibrium tolerance and is independent of the material's TolR.
test NormDispIncr 1.0e-5 20 0
algorithm Newton
pattern Plain 2 Linear {
    load 7 0.0 0.0 -1.0
}

# Retry a failed axial increment as two half increments, following the
# recursive substepping strategy in the published OpenSees triaxial example.
# Each successful substep is committed and recorded; their sum remains exactly
# the requested nominal increment.
set acceptedAxialSubsteps 0
proc advanceDrainedIncrement {increment depth} {
    global maxSubstepDepth acceptedAxialSubsteps
    integrator DisplacementControl 7 3 $increment
    test NormDispIncr 1.0e-5 20 0
    algorithm Newton
    set ok [analyze 1]
    if {$ok == 0} {
        writeState monotonic
        incr acceptedAxialSubsteps
        return 0
    }
    if {$depth >= $maxSubstepDepth} {
        puts stderr "OpenSees drained increment failed at maximum substep depth $depth"
        return -10
    }
    set halfIncrement [expr {0.5*$increment}]
    puts stderr "OpenSees drained increment did not converge; subdividing to $halfIncrement (depth [expr {$depth+1}])"
    set ok [advanceDrainedIncrement $halfIncrement [expr {$depth+1}]]
    if {$ok != 0} {
        return $ok
    }
    return [advanceDrainedIncrement $halfIncrement [expr {$depth+1}]]
}

set nominalIncrement [expr {$finalAxialStrain/$axialSteps}]
for {set i 0} {$i < $axialSteps} {incr i} {
    set ok [advanceDrainedIncrement $nominalIncrement 0]
    if {$ok != 0} {
        puts stderr "OpenSees drained triaxial loading failed at increment $i"
        close $output
        exit 44
    }
}
writeState complete

close $output
puts "OpenSees drained Manzari-Dafalias triaxial history complete: $step states"
puts "Accepted axial increments (including substeps): $acceptedAxialSubsteps"
puts "CSV output: $outputPath"
wipe
exit 0
