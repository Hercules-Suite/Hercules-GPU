# Independent OpenSees characterization of the synthetic Hercules PDMY03
# fixture. Units are Pa, m, and s. These constants are not a calibration.
#
# Usage:
#   OpenSees opensees_pdmy03_cyclic_shear.tcl ?output.csv?

set outputPath "opensees_pdmy03_cyclic_shear.csv"
if {$argc >= 1} { set outputPath [lindex $argv 0] }
if {$argc > 1} {
    puts stderr "usage: OpenSees opensees_pdmy03_cyclic_shear.tcl ?output.csv?"
    exit 2
}

set epsConf [expr {-100000.0/(3.0*2.2e8)}]
set gammaHistory {
     0.0000  0.0002  0.0004  0.0006  0.0008
     0.0006  0.0004  0.0002  0.0000 -0.0002
    -0.0004 -0.0006 -0.0008 -0.0006 -0.0004
    -0.0002  0.0000  0.0000
}

# Domain time 1 establishes 100 kPa isotropic confinement in the forced
# elastic stage. Times 2--19 apply the cyclic shear history after stage 1.
set times {0.0 1.0}
set epsValues [list 0.0 $epsConf]
set gammaValues {0.0 0.0}
set epsGammaValues [list 0.0 $epsConf]
set domainTime 2.0
foreach gamma $gammaHistory {
    lappend times $domainTime
    lappend epsValues $epsConf
    lappend gammaValues $gamma
    lappend epsGammaValues [expr {$epsConf+$gamma}]
    set domainTime [expr {$domainTime+1.0}]
}

wipe
model BasicBuilder -ndm 3 -ndf 3
nDMaterial PressureDependMultiYield03 1 3 0.0 \
    9.0e7 2.2e8 32.0 0.10 1.01e5 0.50 26.0 0 \
    0.067 0.50 0.23 1.0 0.0 0.06 3.0 0.27 \
    20 1.0e4 0.01 1.01e5 1.73e3

node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

# Every displacement is either zero or prescribed by an affine coordinate
# history. ux=epsConf*x+gamma*y, uy=epsConf*y, uz=epsConf*z.
fix 1 1 1 1
fix 2 0 1 1
fix 3 0 0 1
fix 4 0 0 1
fix 5 1 1 0
fix 6 0 1 0
fix 7 0 0 0
fix 8 0 0 0

element stdBrick 1 1 2 3 4 5 6 7 8 1 0.0 0.0 0.0

timeSeries Path 1 -time $times -values $epsValues -useLast
timeSeries Path 2 -time $times -values $gammaValues -useLast
timeSeries Path 3 -time $times -values $epsGammaValues -useLast
pattern Plain 1 1 {
    sp 2 1 1.0
    sp 6 1 1.0
    sp 3 2 1.0
    sp 4 2 1.0
    sp 7 2 1.0
    sp 8 2 1.0
    sp 5 3 1.0
    sp 6 3 1.0
    sp 7 3 1.0
    sp 8 3 1.0
}
pattern Plain 2 2 {
    sp 4 1 1.0
    sp 8 1 1.0
}
pattern Plain 3 3 {
    sp 3 1 1.0
    sp 7 1 1.0
}

constraints Penalty 1.0e18 1.0e18
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-12 50 0
algorithm Newton
integrator LoadControl 1.0
analysis Static

if {[analyze 1] != 0} {
    puts stderr "OpenSees PDMY03 consolidation failed"
    exit 40
}
updateMaterialStage -material 1 -stage 1

set output [open $outputPath w]
puts $output "test_id,driver_version,step,epsilon_xy,gamma_xy,sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_yz,sigma_xz,domain_time"
set step 0
foreach targetGamma $gammaHistory {
    if {[analyze 1] != 0} {
        puts stderr "OpenSees PDMY03 plastic step $step failed"
        close $output
        exit 41
    }
    set stress [eleResponse 1 material 1 stress]
    set strain [eleResponse 1 material 1 strain]
    # PDMY03 3.8 appends pressureD to the six stress components.
    if {[llength $stress] < 6 || [llength $strain] != 6} {
        puts stderr "OpenSees PDMY03 material response has unexpected size at step $step: stress=[llength $stress] strain=[llength $strain]"
        close $output
        exit 42
    }
    set gamma [lindex $strain 3]
    if {abs($gamma-$targetGamma) > 1.0e-12} {
        puts stderr "OpenSees PDMY03 shear path mismatch at step $step"
        close $output
        exit 43
    }
    puts $output [format "pdmy03_cyclic_shear,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e" \
        1 $step [expr {0.5*$gamma}] $gamma \
        [lindex $stress 0] [lindex $stress 1] [lindex $stress 2] \
        [lindex $stress 3] [lindex $stress 4] [lindex $stress 5] [getTime]]
    incr step
}
close $output
puts "OpenSees PDMY03 cyclic shear complete: $step states"
puts "CSV output: $outputPath"
wipe
exit 0
