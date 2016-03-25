source ../misc/KeyFrame3.tcl
source ../misc/eng3d_pillers_textured.tcl

vtkLight newlight
newlight SetPosition 0.75 0.5 1
RenderWindow.renderer AddLight newlight
set lights [RenderWindow.renderer GetLights]
$lights InitTraversal
set first_light [$lights GetNextItem]
$first_light SetIntensity 0

if { [string length [info commands \
	vol]] } {
    volInterface SetCroppingRegionPlanes 0 41.88 0 70 0 72
    volInterface SetCropping 1
}

renderWindowInterface SetCameraFileName ../misc/seq1end.cam
[renderWindow GetRenderer] ResetCameraClippingRange


set number 24
set frame [expr $number - 1]

set p1 [building.p2.actor GetProperty]
set p2 [building.p5.actor GetProperty]

$p1 SetOpacity [expr 1.0 - double($frame)/$number]
$p2 SetOpacity [expr 1.0 - double($frame)/$number]

source ../misc/playmove.tcl
KeyRun playmove 255
