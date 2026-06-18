# OpenSees companion model for the Hercules homogeneous BrickUP column.
#
# Run from this directory:
#   OpenSees model.tcl
#
# The model intentionally keeps the same 3D block geometry as the Hercules
# case instead of imposing equalDOF 1D constraints. Bottom and lateral dashpot
# boundaries are represented with zeroLength viscous dashpots attached to fixed
# anchor nodes.

wipe
model BasicBuilder -ndm 3 -ndf 4

set scriptDir [file dirname [file normalize [info script]]]
set caseDir [file dirname $scriptDir]
set outputDir [file join $scriptDir output]
file mkdir $outputDir

# ----------------------------------------------------------------------
# Geometry and time controls
# ----------------------------------------------------------------------

set lx 1.0
set ly 1.0
set depth 32.0
set nx 1
set ny 1
set nz 32

set dx [expr {$lx / double($nx)}]
set dy [expr {$ly / double($ny)}]
set dz [expr {$depth / double($nz)}]

set dt 0.00005
set tEnd 0.40
set outputEvery 10

# ----------------------------------------------------------------------
# Dry skeleton and pore-fluid properties matching parameters.in
# ----------------------------------------------------------------------

set rho 2000.0
set vs 200.0
set vp 400.0
set fluidBulk 2.2e9
set fluidRho 1000.0
set permX 0.0
set permY 0.0
set permZ 0.0

set shearMod [expr {$rho * $vs * $vs}]
set lambda [expr {$rho * $vp * $vp - 2.0 * $shearMod}]
set young [expr {$shearMod * (3.0 * $lambda + 2.0 * $shearMod) / ($lambda + $shearMod)}]
set poisson [expr {$lambda / (2.0 * ($lambda + $shearMod))}]

puts "Dry skeleton: E=$young Pa, nu=$poisson, rho=$rho kg/m^3"

nDMaterial ElasticIsotropic 1 $young $poisson $rho

# ----------------------------------------------------------------------
# Mesh: 1 x 1 x 32 BrickUP elements. OpenSees z is positive upward, so depth
# is represented as negative z.
# ----------------------------------------------------------------------

proc nodeTag {i j k} {
    return [expr {1 + $i + 2*$j + 4*$k}]
}

array set nodeX {}
array set nodeY {}
array set nodeZ {}

for {set k 0} {$k <= $nz} {incr k} {
    set z [expr {-1.0 * $k * $dz}]
    for {set j 0} {$j <= $ny} {incr j} {
        set y [expr {$j * $dy}]
        for {set i 0} {$i <= $nx} {incr i} {
            set x [expr {$i * $dx}]
            set nd [nodeTag $i $j $k]
            node $nd $x $y $z
            set nodeX($nd) $x
            set nodeY($nd) $y
            set nodeZ($nd) $z
        }
    }
}

set eleTag 1
for {set k 0} {$k < $nz} {incr k} {
    set kb [expr {$k + 1}]
    set n1 [nodeTag 0 0 $kb]
    set n2 [nodeTag 1 0 $kb]
    set n3 [nodeTag 1 1 $kb]
    set n4 [nodeTag 0 1 $kb]
    set n5 [nodeTag 0 0 $k]
    set n6 [nodeTag 1 0 $k]
    set n7 [nodeTag 1 1 $k]
    set n8 [nodeTag 0 1 $k]
    element brickUP $eleTag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 1 \
        $fluidBulk $fluidRho $permX $permY $permZ 0.0 0.0 0.0
    incr eleTag
}

# Drained pore-pressure boundary at the top surface, matching the Hercules
# example setting poroelastic_surface_pressure_bc = drained.
foreach nd [list [nodeTag 0 0 0] [nodeTag 1 0 0] [nodeTag 1 1 0] [nodeTag 0 1 0]] {
    fix $nd 0 0 0 1
}

# ----------------------------------------------------------------------
# Dashpot absorbing boundaries on bottom and lateral faces.
# ----------------------------------------------------------------------

array set dashC {}

proc addDashpotContribution {nd dof c} {
    global dashC
    set key "$nd,$dof"
    if {![info exists dashC($key)]} {
        set dashC($key) 0.0
    }
    set dashC($key) [expr {$dashC($key) + $c}]
}

proc addFaceDashpots {nodes normalDof area rho vp vs} {
    set share [expr {$area / 4.0}]
    foreach nd $nodes {
        for {set dof 1} {$dof <= 3} {incr dof} {
            if {$dof == $normalDof} {
                set speed $vp
            } else {
                set speed $vs
            }
            addDashpotContribution $nd $dof [expr {$rho * $speed * $share}]
        }
    }
}

# Bottom face: normal is vertical.
addFaceDashpots [list [nodeTag 0 0 $nz] [nodeTag 1 0 $nz] \
                      [nodeTag 1 1 $nz] [nodeTag 0 1 $nz]] \
    3 [expr {$dx * $dy}] $rho $vp $vs

# Lateral faces: one quadrilateral face per BrickUP element side.
for {set k 0} {$k < $nz} {incr k} {
    set kt $k
    set kb [expr {$k + 1}]

    addFaceDashpots [list [nodeTag 0 0 $kb] [nodeTag 0 1 $kb] \
                          [nodeTag 0 1 $kt] [nodeTag 0 0 $kt]] \
        1 [expr {$dy * $dz}] $rho $vp $vs
    addFaceDashpots [list [nodeTag 1 0 $kb] [nodeTag 1 1 $kb] \
                          [nodeTag 1 1 $kt] [nodeTag 1 0 $kt]] \
        1 [expr {$dy * $dz}] $rho $vp $vs
    addFaceDashpots [list [nodeTag 0 0 $kb] [nodeTag 1 0 $kb] \
                          [nodeTag 1 0 $kt] [nodeTag 0 0 $kt]] \
        2 [expr {$dx * $dz}] $rho $vp $vs
    addFaceDashpots [list [nodeTag 0 1 $kb] [nodeTag 1 1 $kb] \
                          [nodeTag 1 1 $kt] [nodeTag 0 1 $kt]] \
        2 [expr {$dx * $dz}] $rho $vp $vs
}

set nextNodeTag 100000
set nextMatTag 100000
set nextEleTag 100000

foreach key [array names dashC] {
    set c $dashC($key)
    if {$c <= 0.0} {
        continue
    }
    set parts [split $key ","]
    set nd [lindex $parts 0]
    set dof [lindex $parts 1]

    set anchor $nextNodeTag
    incr nextNodeTag
    node $anchor $nodeX($nd) $nodeY($nd) $nodeZ($nd)
    fix $anchor 1 1 1 1

    set mat $nextMatTag
    incr nextMatTag
    uniaxialMaterial Viscous $mat $c 1.0

    element zeroLength $nextEleTag $nd $anchor -mat $mat -dir $dof
    incr nextEleTag
}

puts "Created [array size dashC] dashpot coefficients on bottom/lateral boundaries"

# ----------------------------------------------------------------------
# Base force. The same velocity history used by Hercules is converted to
# total base force by rho_base * Vs_base * A_base * velocity(t), then split by
# bottom-node tributary area.
# ----------------------------------------------------------------------

proc readScalarSeries {path} {
    set fp [open $path r]
    set n [string trim [gets $fp]]
    set t0 [string trim [gets $fp]]
    set dtSeries [string trim [gets $fp]]
    set values {}
    for {set i 0} {$i < $n} {incr i} {
        if {[gets $fp line] < 0} {
            close $fp
            error "Unexpected end of $path while reading sample $i of $n"
        }
        lappend values [string trim $line]
    }
    close $fp
    return [list $n $t0 $dtSeries $values]
}

set seriesPath [file join $caseDir sourcefiles baseforcefunction.in]
lassign [readScalarSeries $seriesPath] nSeries t0Series dtSeries velocitySeries
if {[expr {abs($t0Series)}] > 1.0e-12} {
    puts "WARNING: nonzero base-force t0=$t0Series is not shifted in this Tcl model"
}

set rhoBase 2000.0
set vsBase 200.0
set baseArea [expr {$lx * $ly}]
set baseForceScale 1.0
set totalForceFactor [expr {$baseForceScale * $rhoBase * $vsBase * $baseArea}]

set dirX 1.0
set dirY 0.0
set dirZ 0.0

timeSeries Path 1 -dt $dtSeries -values $velocitySeries -factor $totalForceFactor

set bottomNodes [list [nodeTag 0 0 $nz] [nodeTag 1 0 $nz] \
                      [nodeTag 1 1 $nz] [nodeTag 0 1 $nz]]
pattern Plain 1 1 {
    foreach nd $bottomNodes {
        load $nd [expr {0.25 * $dirX}] [expr {0.25 * $dirY}] [expr {0.25 * $dirZ}] 0.0
    }
}

# ----------------------------------------------------------------------
# Station output at the same depths as the Hercules example. Values are
# trilinear averages at x = 0.5 m, y = 0.5 m.
# ----------------------------------------------------------------------

set stationDepths [list 0.0 8.0 16.0 24.0 31.5]
array set stationFile {}

for {set s 0} {$s < [llength $stationDepths]} {incr s} {
    set stationFile($s) [open [file join $outputDir "station.$s"] w]
    puts $stationFile($s) "#  Time(s)         X|(m)         Y-(m)         Z.(m)       X|(m/s)       Y-(m/s)       Z.(m/s)      X|(m/s2)      Y-(m/s2)      Z.(m/s2)      Pressure(Pa)"
}

proc planeAverage {k dof quantity} {
    set nodes [list [nodeTag 0 0 $k] [nodeTag 1 0 $k] [nodeTag 1 1 $k] [nodeTag 0 1 $k]]
    set sum 0.0
    foreach nd $nodes {
        if {$quantity eq "disp"} {
            set value [nodeDisp $nd $dof]
        } elseif {$quantity eq "vel"} {
            set value [nodeVel $nd $dof]
        } elseif {$quantity eq "accel"} {
            set value [nodeAccel $nd $dof]
        } else {
            error "Unknown station quantity: $quantity"
        }
        set sum [expr {$sum + $value}]
    }
    return [expr {$sum / 4.0}]
}

proc depthAverage {depth dof quantity} {
    global dz nz
    set scaled [expr {$depth / $dz}]
    set k0 [expr {int(floor($scaled))}]
    if {$k0 < 0} {
        set k0 0
        set xi 0.0
    } elseif {$k0 >= $nz} {
        set k0 [expr {$nz - 1}]
        set xi 1.0
    } else {
        set xi [expr {$scaled - $k0}]
    }
    set k1 [expr {$k0 + 1}]
    set v0 [planeAverage $k0 $dof $quantity]
    set v1 [planeAverage $k1 $dof $quantity]
    return [expr {(1.0 - $xi) * $v0 + $xi * $v1}]
}

proc writeStations {time} {
    global stationDepths stationFile
    for {set s 0} {$s < [llength $stationDepths]} {incr s} {
        set depth [lindex $stationDepths $s]
        set ux [depthAverage $depth 1 disp]
        set uy [depthAverage $depth 2 disp]
        set uz [depthAverage $depth 3 disp]
        set vx [depthAverage $depth 1 vel]
        set vy [depthAverage $depth 2 vel]
        set vz [depthAverage $depth 3 vel]
        set ax [depthAverage $depth 1 accel]
        set ay [depthAverage $depth 2 accel]
        set az [depthAverage $depth 3 accel]
        set pp [depthAverage $depth 4 disp]
        puts $stationFile($s) [format "%10.6f % .6e % .6e % .6e % .6e % .6e % .6e % .6e % .6e % .6e % .6e" \
            $time $ux $uy $uz $vx $vy $vz $ax $ay $az $pp]
        flush $stationFile($s)
    }
}

# ----------------------------------------------------------------------
# Transient analysis
# ----------------------------------------------------------------------

constraints Transformation
numberer RCM
if {[catch {system UmfPack} solverError]} {
    puts "UmfPack unavailable ($solverError); using SparseGeneral -piv"
    system SparseGeneral -piv
}
test NormDispIncr 1.0e-10 20 0
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set nSteps [expr {int(round($tEnd / $dt))}]
writeStations 0.0

set ok 0
for {set step 1} {$step <= $nSteps} {incr step} {
    set ok [analyze 1 $dt]
    if {$ok != 0} {
        puts "OpenSees analysis failed at step $step time [getTime] with code $ok"
        break
    }
    if {[expr {$step % $outputEvery}] == 0} {
        writeStations [getTime]
    }
}

for {set s 0} {$s < [llength $stationDepths]} {incr s} {
    close $stationFile($s)
}

if {$ok == 0} {
    puts "OpenSees BrickUP column complete: outputs in $outputDir"
}
