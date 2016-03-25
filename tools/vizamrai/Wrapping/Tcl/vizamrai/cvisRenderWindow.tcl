##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisRenderWindow.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Render Window for 3D data
##

proc cvisTkCon {} {
    set tracervar [::itcl::local cvisTrace #auto]
    if {[info exists ::tkcon::PRIV(root)]} {
	if [winfo ismapped $::tkcon::PRIV(root)] {
	    tkcon hide
	} {
	    tkcon show
	}
    } {
	global argv argc
	# tkCon is reading command line args, this is bad
	set argv ""
	set argc 0

	font create vizamraitkconn -family Courier -size 10
	tkcon font vizamraitkconn
	tkcon attach Main
    }
}

class cvisView {
    variable Views
    variable ViewNames

    constructor {} {
	::cvis::RestartFile::Add $this
    }

    method SetClippingRange {name min max} {
	set ViewNames($name) 1
	set Views($name,ClippingRange) "$min $max"
    }
    
    method SetDistance {name distance} {
	set ViewNames($name) 1
	set Views($name,Distance) $distance
    }

    method SetEyeAngle {name angle} {
	set ViewNames($name) 1
	set Views($name,EyeAngle) $angle
    }

    method SetFocalDisk {name focaldisk} {
	set ViewNames($name) 1
	set Views($name,FocalDisk) $focaldisk
    }


    method SetFocalPoint {name x y z} {
	set ViewNames($name) 1
	set Views($name,FocalPoint) "$x $y $z"
    }

    method SetPosition {name x y z} {
	set ViewNames($name) 1
	set Views($name,Position) "$x $y $z"
    }

    method SetParallelProjection {name parallelprojection} {
	set ViewNames($name) 1
	set Views($name,ParallelProjection) $parallelprojection
    }

    method SetThickness {name thickness} {
	set Views($name,Thickness) $thickness
    }

    method SetViewAngle {name viewangle} {
	set Views($name,ViewAngle) $viewangle
    }

    method SetViewUp {name x y z} {
	set Views($name,ViewUp) "$x $y $z"
    }

    method Save {file} {
	foreach name [array names ViewNames] {
	    foreach i "ClippingRange Distance EyeAngle FocalDisk FocalPoint ParallelProjection Position Thickness ViewAngle ViewUp" {
		puts $file "$this Set$i $name $Views($name,$i)"
	    }
	}
    }

    # Sets view based on a camera
    method SetViewFromCamera {name camera} {
	set ViewNames($name) 1
	set Views($name,ClippingRange) [$camera GetClippingRange]
	set Views($name,Distance)      [$camera GetDistance]
	set Views($name,EyeAngle)      [$camera GetEyeAngle]
	set Views($name,FocalDisk)     [$camera GetFocalDisk]
	set Views($name,FocalPoint)    [$camera GetFocalPoint]
	set Views($name,Position)      [$camera GetPosition]
	set Views($name,ParallelProjection) \
		[expr int([$camera GetParallelProjection])]
	set Views($name,Thickness)     [$camera GetThickness]
	set Views($name,ViewAngle)     [$camera GetViewAngle]
	set Views($name,ViewUp)        [$camera GetViewUp]
    }

    method GetViewToCamera {name camera} {
	if {[array exists ViewNames]} {
	    if {[array names ViewNames $name] > 0 } {
	    eval $camera SetClippingRange $Views($name,ClippingRange)
	    $camera SetDistance	 $Views($name,Distance)
	    $camera SetEyeAngle	 $Views($name,EyeAngle) 
	    $camera SetFocalDisk	 $Views($name,FocalDisk)
	    eval $camera SetFocalPoint	 $Views($name,FocalPoint)
	    eval $camera SetPosition	 $Views($name,Position)
	    $camera SetParallelProjection	$Views($name,ParallelProjection) 
	    $camera SetThickness	 $Views($name,Thickness)  
	    $camera SetViewAngle	 $Views($name,ViewAngle)   
	    eval $camera SetViewUp	 $Views($name,ViewUp)
	    return 1
	    }
	}
	return 0
    }

    method GetViewNames {} {
	return [array names ViewNames]
    }
}


class cvisRenderWindow {
    variable Name
    variable Create 1
    
    variable UniformScaling 0

    variable ParallelProjection 0
    
    variable WindowNameProp
    
    variable Renderer ""
    variable RenderWindow ""
    
    variable RenderWidget ""

    variable ImageWriter ""
    
    variable Size
    
    # Interactors
    variable InExpose 0
    
    variable RubberZoomPerformed 0
    
    variable RendererFound 0
    
    variable LastX 0
    variable LastY 0
    
    variable WindowCenterX 0 
    variable WindowCenterY 0
    
    variable TkInteractor_InteractiveUpdateRate 15.0
    variable TkInteractor_StillUpdateRate 1.0E-37

    variable InteractiveUpdateRate 15.0
    variable StillUpdateRate 1.0E-37

    # Rubber band zoom objects

    variable RubberBandPoints   ""
    variable RubberBandLines    ""
    variable RubberBandScalars  ""
    variable RubberBandPolyData ""
    variable RubberBandMapper   ""
    variable RubberBandActor    ""
    variable RubberBandColors   ""

    # State of the rubber band zoom
    variable StartRubberZoomX 0
    variable StartRubberZoomY 0
    variable EndRubberZoomX 0
    variable EndRubberZoomY 0

    variable ComputeBoundsFlag 1
    variable Bounds 

    variable Interface ""

    method SetInterface {interface} {
	set Interface $interface
    }

    method GetInterface {} {
	return $Interface
    }

    method Save {file} {
	set camera [$Renderer GetActiveCamera]
	foreach i "ClippingRange Distance EyeAngle FocalDisk FocalPoint ParallelProjection Position Thickness ViewAngle ViewUp" {
	    puts $file "$this Set$i [$camera Get$i]"
	}
    }

    method SetClippingRange {min max} {
	[$Renderer GetActiveCamera] SetClippingRange $min $max
    }

    method SetDistance { distance} {
	[$Renderer GetActiveCamera] SetDistance $distance
    }

    method SetEyeAngle { angle} {
	[$Renderer GetActiveCamera] SetEyeAngle $angle
    }

    method SetFocalDisk { focaldisk} {
	[$Renderer GetActiveCamera] SetFocalDisk $focaldisk
    }

    method SetFocalPoint { x y z} {
	[$Renderer GetActiveCamera] SetFocalPoint $x $y $z
    }

    method SetPosition { x y z} {
	[$Renderer GetActiveCamera] SetPosition $x $y $z
    }

    method SetParallelProjection { parallelprojection} {
	[$Renderer GetActiveCamera] SetParallelProjection $parallelprojection
    }

    method SetThickness { thickness} {
	[$Renderer GetActiveCamera] SetThickness $thickness
    }

    method SetViewAngle { viewangle} {
	[$Renderer GetActiveCamera] SetViewAngle $viewangle
    }

    method SetViewUp { x y z} {
	[$Renderer GetActiveCamera] SetViewUp $x $y $z
    }

    method UpdateRenderer {x y} {
	set tracervar [::itcl::local cvisTrace #auto]

	update
	
	# Get the renderer window dimensions
	set WindowX [lindex [$RenderWidget configure -width] 4]
	set WindowY [lindex [$RenderWidget configure -height] 4]
	
	# Find which renderer event has occurred in
	set renderers [$RenderWindow GetRenderers]
	set numRenderers [$renderers GetNumberOfItems]

	# SGS removed multi renderer capablity
	$renderers InitTraversal; 
	set CurrentRenderer [$renderers GetNextItem]
	set viewport [$CurrentRenderer GetViewport]
	set vpxmin [lindex $viewport 0]
	set vpymin [lindex $viewport 1]
	set vpxmax [lindex $viewport 2]
	set vpymax [lindex $viewport 3]
	set RendererFound 1
	set WindowCenterX [expr double($WindowX)*(($vpxmax + $vpxmin)/2.0)]
	set WindowCenterY [expr double($WindowY)*(($vpymax + $vpymin)/2.0)]
	
	set CurrentCamera [$Renderer GetActiveCamera]
	
	set LastX $x
	set LastY $y
    }
    
    variable CurrentCamera ""
    
    method StartPan {widget x y} {
	$this StartMotion $widget $x $y
	$widget select
    }

    method EndPan {widget x y} {
	$this EndMotion $widget $x $y
	$widget deselect
    }

    method StartRotate {widget x y} {
	$this StartMotion $widget $x $y
	$widget select
    }

    method EndRotate {widget x y} {
	$this EndMotion $widget $x $y
	$widget deselect
    }

    method StartZoom {widget x y} {
	$this StartMotion $widget $x $y
	$widget select
    }

    method EndZoom {widget x y} {
	$this EndMotion $widget $x $y
	$widget deselect
    }
    
    method StartMotion {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	UpdateRenderer $x $y
	
	if { $RendererFound } { 
	    set RubberZoomPerformed 0
	    
	    $RenderWindow SetDesiredUpdateRate $InteractiveUpdateRate
	}
    }
    
    method EndMotion {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]
	if { $RendererFound } {
	    $RenderWindow SetDesiredUpdateRate $StillUpdateRate
	    
	    if { $RubberZoomPerformed } {
		$Renderer RemoveProp $RubberBandActor
		DoRubberZoom $widget
	    }
	    
	    RenderInteractive
	}
    }
    
    method Wireframe {} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	set actors [$Renderer GetActors]
	
	$actors InitTraversal
	set actor [$actors GetNextItem]
	while { $actor != "" } {
	    [$actor GetProperty] SetRepresentationToWireframe
	    set actor [$actors GetNextItem]
	}

	Render
    }
    
    method Surface {} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	set actors [$Renderer GetActors]
	
	$actors InitTraversal
	set actor [$actors GetNextItem]
	while { $actor != "" } {
	    [$actor GetProperty] SetRepresentationToSurface
	    set actor [$actors GetNextItem]
	}

	Render
    }
    
    variable PickedAssembly ""
    variable PrePickedProperty ""

    method PickActor {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	set WindowY [lindex [$widget configure -height] 4]
	
	if { ! $RendererFound } { return }

	if { $PickedAssembly != "" && $PrePickedProperty != "" } {
	    $PickedAssembly SetProperty $PrePickedProperty
	    # release hold on the property
	    $PrePickedProperty UnRegister $PrePickedProperty
	    set PrePickedProperty ""
	}

	if { $assembly != "" } {
	    set PickedAssembly $assembly
	    set PrePickedProperty [$PickedAssembly GetProperty]
	    # hold onto the property
	    $PrePickedProperty Register $PrePickedProperty
	}
	
	RenderInteractive
    }

    method Reset {} {
	$Renderer ResetCamera
	RenderInteractive
    }

    method ResetInteractorFullScreen {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]
	vtkMath math
	set PI [math Pi]
	
	ComputeBounds
	
	set camera [$Renderer GetActiveCamera]
	
	set vn(x) [expr double([lindex [$camera GetViewPlaneNormal] 0])]
	set vn(y) [expr double([lindex [$camera GetViewPlaneNormal] 1])]
	set vn(z) [expr double([lindex [$camera GetViewPlaneNormal] 2])]
	
  	set center(x) [expr ($Bounds(min,x) + $Bounds(max,x))/2.0]
  	set center(y) [expr ($Bounds(min,y) + $Bounds(max,y))/2.0]
  	set center(z) [expr ($Bounds(min,z) + $Bounds(max,z))/2.0]

  	set pscale [expr ($Bounds(max,x) - $Bounds(min,x)) / 2.0 ]
  	if { $pscale > [expr ($Bounds(max,y) - $Bounds(min,y)) / 2.0] } {
  	    set pscale [expr ($Bounds(max,y) - $Bounds(min,y)) / 2.0]
  	}
	
  	set width [expr $Bounds(max,y) - $Bounds(min,y)]
  	if { $width < [expr $Bounds(max,x) - $Bounds(min,x)] } {
  	    set width [expr $Bounds(max,x) - $Bounds(min,x)]
  	}

	set viewangle [expr double([$camera GetViewAngle])]
  	set distance [expr 0.8 * $width / tan($viewangle * $PI/360.0)]

  	set distance [expr $distance + ($Bounds(max,z) - $Bounds(min,z))/2.0]

	set vup(x) [expr double([lindex [$camera GetViewUp] 0])]
	set vup(y) [expr double([lindex [$camera GetViewUp] 1])]
	set vup(z) [expr double([lindex [$camera GetViewUp] 2])]
	
  	set dot [expr ($vup(x) * $vn(x)) + ($vup(y) * $vn(y)) +  ($vup(z) * $vn(z))]
  	if {$dot > 0.999} {
  	    $camera SetViewUp [expr -$vup(z)] $vup(0) $vup(1)
  	}
	
  	$camera SetFocalPoint $center(x) $center(y) $center(z)
  	$camera SetPosition \
  		[expr $center(x) + $distance*$vn(x)] \
  		[expr $center(y) + $distance*$vn(y)] \
  		[expr $center(z) + $distance*$vn(z)]
  	$Renderer ResetCameraClippingRange

  	Render

	set center(x) [expr double([lindex [$camera GetFocalPoint] 0])]
	set center(y) [expr double([lindex [$camera GetFocalPoint] 1])]
	set center(z) [expr double([lindex [$camera GetFocalPoint] 2])]

	set position(x) [expr double([lindex [$camera GetPosition] 0])]
	set position(y) [expr double([lindex [$camera GetPosition] 1])]
	set position(z) [expr double([lindex [$camera GetPosition] 2])]
	
	$Renderer SetWorldPoint $center(x) $center(y) $center(z) 1.0
	$Renderer WorldToDisplay
	set center_screen(x) [lindex [$Renderer GetDisplayPoint] 0]
	set center_screen(y) [lindex [$Renderer GetDisplayPoint] 1]
	set center_screen(z) [lindex [$Renderer GetDisplayPoint] 2]

	set old_pscale [$camera GetParallelScale]
	
	foreach plane "x y" {
	    set Screen(max,$plane) [expr -$::cvis::math::FLT_MAX]
	}
	
	set screen_size(x) [expr double([lindex [$Renderer GetSize] 0])]
	set screen_size(y) [expr double([lindex [$Renderer GetSize] 1])]

	foreach bounds_x {min max min max min max min max} bounds_y {min min max max min min max max} bounds_z {min min min min max max max max} {
	    
	    $Renderer SetWorldPoint $Bounds($bounds_x,x) $Bounds($bounds_y,y)\
		    $Bounds($bounds_z,z) 1.0
	    $Renderer WorldToDisplay
	    set x [lindex [$Renderer GetDisplayPoint] 0]
	    set y [lindex [$Renderer GetDisplayPoint] 1]
	    set z [lindex [$Renderer GetDisplayPoint] 2]
	    
	    if { $Screen(max,x) < abs($x - $screen_size(x) / 2.0) } {
		set Screen(max,x) [expr abs($x - $screen_size(x) / 2.0)]
		set Screen(max,zx) $z
	    }
	    
	    if { $Screen(max,y) < abs($y - $screen_size(y) / 2.0) } {
		set Screen(max,y) [expr abs($y - $screen_size(y) / 2.0)]
		set Screen(max,zy) $z
	    }
	}

	set pscale [expr double(2*$Screen(max,x)) / $screen_size(x)]
	if { $pscale < [expr double(2*$Screen(max,y)) / $screen_size(y)] } {
	    set pscale [expr double(2*$Screen(max,y)) / $screen_size(y)]

	    $Renderer SetDisplayPoint \
		    [expr $screen_size(x)/2.0] \
		    [expr $screen_size(y)/2.0] \
		    $Screen(max,zy)
	    $Renderer DisplayToWorld
	    set nearFocalPoint [$Renderer GetWorldPoint]
	    set nearFocalPoint0 [lindex $nearFocalPoint 0]
	    set nearFocalPoint1 [lindex $nearFocalPoint 1]
	    set nearFocalPoint2 [lindex $nearFocalPoint 2]
	    set nearFocalPoint3 [lindex $nearFocalPoint 3]
	    if { $nearFocalPoint3 != 0.0 } {
	    set nearFocalPoint0 [expr $nearFocalPoint0 / $nearFocalPoint3]
		set nearFocalPoint1 [expr $nearFocalPoint1 / $nearFocalPoint3]
		set nearFocalPoint2 [expr $nearFocalPoint2 / $nearFocalPoint3]
	    }

	    $Renderer SetDisplayPoint \
		    [expr $screen_size(x)/2.0] \
		    [expr $screen_size(y)/2.0 + $Screen(max,y)] \
		    $Screen(max,zy)
	    $Renderer DisplayToWorld
	    set nearFocalEdge [$Renderer GetWorldPoint]
	    set nearFocalEdge0 [expr double([lindex $nearFocalEdge 0])]
	    set nearFocalEdge1 [expr double([lindex $nearFocalEdge 1])]
	    set nearFocalEdge2 [expr double([lindex $nearFocalEdge 2])]
	    set nearFocalEdge3 [expr double([lindex $nearFocalEdge 3])]
	    if { $nearFocalEdge3 != 0.0 } {
		set nearFocalEdge0 [expr $nearFocalEdge0 / $nearFocalEdge3]
		set nearFocalEdge1 [expr $nearFocalEdge1 / $nearFocalEdge3]
		set nearFocalEdge2 [expr $nearFocalEdge2 / $nearFocalEdge3]
	    }

	    $Renderer SetDisplayPoint \
		    [expr $screen_size(x)/2.0] \
		    0.0 \
		    $center_screen(z)
	    $Renderer DisplayToWorld
	    set projectionPlaneEdge [$Renderer GetWorldPoint]
	    set projectionPlaneEdge0 [expr double([lindex $projectionPlaneEdge 0])]
	    set projectionPlaneEdge1 [expr double([lindex $projectionPlaneEdge 1])]
	    set projectionPlaneEdge2 [expr double([lindex $projectionPlaneEdge 2])]
	    set projectionPlaneEdge3 [expr double([lindex $projectionPlaneEdge 3])]
	    if { $projectionPlaneEdge3 != 0.0 } {
		set projectionPlaneEdge0 [expr $projectionPlaneEdge0 / $projectionPlaneEdge3]
		set projectionPlaneEdge1 [expr $projectionPlaneEdge1 / $projectionPlaneEdge3]
		set projectionPlaneEdge2 [expr $projectionPlaneEdge2 / $projectionPlaneEdge3]
	    }
	} {
	    $Renderer SetDisplayPoint \
		    [expr $screen_size(x)/2.0] \
		    [expr $screen_size(y)/2.0] \
		    $Screen(max,zx)
	    $Renderer DisplayToWorld
	    set nearFocalPoint [$Renderer GetWorldPoint]
	    set nearFocalPoint0 [lindex $nearFocalPoint 0]
	    set nearFocalPoint1 [lindex $nearFocalPoint 1]
	    set nearFocalPoint2 [lindex $nearFocalPoint 2]
	    set nearFocalPoint3 [lindex $nearFocalPoint 3]
	    if { $nearFocalPoint3 != 0.0 } {
	    set nearFocalPoint0 [expr $nearFocalPoint0 / $nearFocalPoint3]
		set nearFocalPoint1 [expr $nearFocalPoint1 / $nearFocalPoint3]
		set nearFocalPoint2 [expr $nearFocalPoint2 / $nearFocalPoint3]
	    }

	    $Renderer SetDisplayPoint \
		    [expr $screen_size(x)/2.0 + $Screen(max,x)] \
		    [expr $screen_size(y)/2.0] \
		    $Screen(max,zx)
	    $Renderer DisplayToWorld
	    set nearFocalEdge [$Renderer GetWorldPoint]
	    set nearFocalEdge0 [expr double([lindex $nearFocalEdge 0])]
	    set nearFocalEdge1 [expr double([lindex $nearFocalEdge 1])]
	    set nearFocalEdge2 [expr double([lindex $nearFocalEdge 2])]
	    set nearFocalEdge3 [expr double([lindex $nearFocalEdge 3])]
	    if { $nearFocalEdge3 != 0.0 } {
		set nearFocalEdge0 [expr $nearFocalEdge0 / $nearFocalEdge3]
		set nearFocalEdge1 [expr $nearFocalEdge1 / $nearFocalEdge3]
		set nearFocalEdge2 [expr $nearFocalEdge2 / $nearFocalEdge3]
	    }


	    $Renderer SetDisplayPoint \
		    0.0 \
		    [expr $screen_size(y)/2.0] \
		    $center_screen(z)
	    $Renderer DisplayToWorld
	    set projectionPlaneEdge [$Renderer GetWorldPoint]
	    set projectionPlaneEdge0 [expr double([lindex $projectionPlaneEdge 0])]
	    set projectionPlaneEdge1 [expr double([lindex $projectionPlaneEdge 1])]
	    set projectionPlaneEdge2 [expr double([lindex $projectionPlaneEdge 2])]
	    set projectionPlaneEdge3 [expr double([lindex $projectionPlaneEdge 3])]
	    if { $projectionPlaneEdge3 != 0.0 } {
		set projectionPlaneEdge0 [expr $projectionPlaneEdge0 / $projectionPlaneEdge3]
		set projectionPlaneEdge1 [expr $projectionPlaneEdge1 / $projectionPlaneEdge3]
		set projectionPlaneEdge2 [expr $projectionPlaneEdge2 / $projectionPlaneEdge3]
	    }
	}

	
	#
	set opp_dist [expr \
		sqrt( \
		($projectionPlaneEdge0 - $center(x))*($projectionPlaneEdge0 - $center(x)) + \
		($projectionPlaneEdge1 - $center(y))*($projectionPlaneEdge1 - $center(y)) + \
		($projectionPlaneEdge2 - $center(z))*($projectionPlaneEdge2 - $center(z)) )]

	set adj_dist [expr \
		sqrt( \
		($position(x) - $center(x))*($position(x) - $center(x)) + \
		($position(y) - $center(y))*($position(y) - $center(y)) + \
		($position(z) - $center(z))*($position(z) - $center(z)) )]

	set angle [expr atan($opp_dist / $adj_dist)]

	set ydist [expr \
		sqrt( \
		($nearFocalPoint0 - $nearFocalEdge0)*($nearFocalPoint0 - $nearFocalEdge0) + \
		($nearFocalPoint1 - $nearFocalEdge1)*($nearFocalPoint1 - $nearFocalEdge1) + \
		($nearFocalPoint2 - $nearFocalEdge2)*($nearFocalPoint2 - $nearFocalEdge2) )]

	set d [expr $ydist/tan($angle)]

	set nearDist [expr \
		sqrt( \
		($nearFocalPoint0 - $center(x))*($nearFocalPoint0 - $center(x)) + \
		($nearFocalPoint1 - $center(y))*($nearFocalPoint1 - $center(y)) + \
		($nearFocalPoint2 - $center(z))*($nearFocalPoint2 - $center(z)) )]

	set factor [expr [$camera GetDistance] / ($d + $nearDist)]
	$CurrentCamera Dolly $factor	

	# I do this better than rubber band zoom!
	$camera SetParallelScale [expr $old_pscale * $pscale]

  	$Renderer ResetCameraClippingRange
	Render

	math Delete
    }

    method ResetInteractor {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]

	# Get the renderer window dimensions
	set WindowX [lindex [$widget configure -width] 4]
	set WindowY [lindex [$widget configure -height] 4]

	# Find which renderer event has occurred in
	set CurrentRenderWindow [$widget GetRenderWindow]
	set renderers [$CurrentRenderWindow GetRenderers]
	set numRenderers [$renderers GetNumberOfItems]

	$renderers InitTraversal; set RendererFound 0
	for {set i 0} {$i < $numRenderers} {incr i} {
	    set CurrentRenderer [$renderers GetNextItem]
	    set vx [expr double($x) / $WindowX]
	    set vy [expr ($WindowY - double($y)) / $WindowY]

	    set viewport [$CurrentRenderer GetViewport]
	    set vpxmin [lindex $viewport 0]
	    set vpymin [lindex $viewport 1]
	    set vpxmax [lindex $viewport 2]
	    set vpymax [lindex $viewport 3]
	    if { $vx >= $vpxmin && $vx <= $vpxmax && \
		    $vy >= $vpymin && $vy <= $vpymax} {
		set RendererFound 1
		break
	    }
	}

	if { $RendererFound } {$CurrentRenderer ResetCamera}

	Render
    }

    method HideSideToolbar {widget} {
	if {[winfo ismapped $widget.tools]} {
	    pack forget $widget.tools
	    pack $widget.empty
	} {
	    pack forget $widget.empty 
	    pack $widget.tools 
	}
    }

    method HideBottomToolbar {widget} {
	if {[winfo ismapped $widget.tools]} {
	    pack forget $widget.tools
	    pack $widget.empty 
	} {
	    pack forget $widget.empty 
	    pack $widget.tools -fill x 
	}
    }

    method Rotate {x y} {
	set tracervar [::itcl::local cvisTrace #auto]

	if { $RendererFound } {
	    $CurrentCamera Azimuth [expr ($LastX - $x)]
	    $CurrentCamera Elevation [expr ($y - $LastY)]
	    $CurrentCamera OrthogonalizeViewUp
	    $Renderer ResetCameraClippingRange

	    set LastX $x
	    set LastY $y

	    RenderInteractive
	}
    }

    method Pan {x y} {
	set tracervar [::itcl::local cvisTrace #auto]

	if { ! $RendererFound } { return }

	set FPoint [$CurrentCamera GetFocalPoint]
        set FPoint0 [lindex $FPoint 0]
        set FPoint1 [lindex $FPoint 1]
        set FPoint2 [lindex $FPoint 2]
	
	set PPoint [$CurrentCamera GetPosition]
        set PPoint0 [lindex $PPoint 0]
        set PPoint1 [lindex $PPoint 1]
        set PPoint2 [lindex $PPoint 2]
	
	$Renderer SetWorldPoint $FPoint0 $FPoint1 $FPoint2 1.0
	$Renderer WorldToDisplay
	set DPoint [$Renderer GetDisplayPoint]
	set focalDepth [lindex $DPoint 2]

	set APoint0 [expr $WindowCenterX + ($x - $LastX)]
	set APoint1 [expr $WindowCenterY - ($y - $LastY)]
	
	$Renderer SetDisplayPoint $APoint0 $APoint1 $focalDepth
	$Renderer DisplayToWorld
	set RPoint [$Renderer GetWorldPoint]
        set RPoint0 [lindex $RPoint 0]
        set RPoint1 [lindex $RPoint 1]
        set RPoint2 [lindex $RPoint 2]
        set RPoint3 [lindex $RPoint 3]
	if { $RPoint3 != 0.0 } {
	    set RPoint0 [expr $RPoint0 / $RPoint3]
	    set RPoint1 [expr $RPoint1 / $RPoint3]
	    set RPoint2 [expr $RPoint2 / $RPoint3]
	}
	
	$CurrentCamera SetFocalPoint \
		[expr ($FPoint0 - $RPoint0)/2.0 + $FPoint0] \
		[expr ($FPoint1 - $RPoint1)/2.0 + $FPoint1] \
		[expr ($FPoint2 - $RPoint2)/2.0 + $FPoint2]
	
	$CurrentCamera SetPosition \
		[expr ($FPoint0 - $RPoint0)/2.0 + $PPoint0] \
		[expr ($FPoint1 - $RPoint1)/2.0 + $PPoint1] \
		[expr ($FPoint2 - $RPoint2)/2.0 + $PPoint2]

	set LastX $x
	set LastY $y
	
	RenderInteractive
    }

    method WheelZoom {d} {
	set tracervar [::itcl::local cvisTrace #auto]

	if { ! $RendererFound } {
	    return 
	}

	set zoomFactor [expr pow(1.02,$d/20)]
	
	$this DoZoom $zoomFactor
    }

    method Zoom {x y} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	if { ! $RendererFound } {
	    return 
	}
	
	set zoomFactor [expr pow(1.02,(0.5*($y - $LastY)))]
	
	$this DoZoom $zoomFactor
	
	set LastX $x
	set LastY $y
    }

    method DoZoom {zoomFactor} {
	if {[$CurrentCamera GetParallelProjection]} {
	    set parallelScale [expr [$CurrentCamera GetParallelScale] * $zoomFactor];
	    $CurrentCamera SetParallelScale $parallelScale;
	} else {
	    $CurrentCamera Dolly $zoomFactor
	    $Renderer ResetCameraClippingRange
	}

	RenderInteractive
    }

    method Top {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) $FPoint(x)
	set CPoint(y) $FPoint(y)
	set CPoint(z) [expr $FPoint(z) + 10.0]

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 1 0

	Reset
    }

    method Bottom {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) $FPoint(x)
	set CPoint(y) $FPoint(y)
	set CPoint(z) [expr $FPoint(z) - 10.0]

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 -1 0

	Reset
    }

    method Left {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) [expr $FPoint(x) - 10.0]
	set CPoint(y) $FPoint(y)
	set CPoint(z) $FPoint(z)

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 0 1

	Reset
    }

    method Right {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) [expr $FPoint(x) + 10.0]
	set CPoint(y) $FPoint(y)
	set CPoint(z) $FPoint(z)

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 0 1

	Reset
    }

    method Front {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) $FPoint(x)
	set CPoint(y) [expr $FPoint(y) + 10.0]
	set CPoint(z) $FPoint(z)

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 0 1

	Reset
    }

    method Back {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) $FPoint(x)
	set CPoint(y) [expr $FPoint(y) - 10.0]
	set CPoint(z) $FPoint(z)

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 0 1

	Reset
    }
    
    method BirdsEye {} {
	ComputeBounds

	set FPoint(x) [expr ($Bounds(max,x) + $Bounds(min,x)) / 2.0 ]
	set FPoint(y) [expr ($Bounds(max,y) + $Bounds(min,y)) / 2.0 ]
	set FPoint(z) [expr ($Bounds(max,z) + $Bounds(min,z)) / 2.0 ]
	
	set CPoint(x) [expr $FPoint(x) + 10.0]
	set CPoint(y) [expr $FPoint(y) + 10.0]
	set CPoint(z) [expr $FPoint(z) + 10.0]

	set camera [$Renderer GetActiveCamera]
	$camera SetFocalPoint $FPoint(x) $FPoint(y) $FPoint(z) 
	$camera SetPosition   $CPoint(x) $CPoint(y) $CPoint(z) 
	$camera SetViewUp 0 0 1

	Reset
    }

    method SetView {name} {
	$this.views SetViewFromCamera $name [$Renderer GetActiveCamera]
    }

    method GetView {name} {
	# Get View if it exists otherwise default to birds eye
	if {[$this.views GetViewToCamera $name [$Renderer GetActiveCamera]]} {
	    Render
	} {
	    BirdsEye
	}

    }

    method GetBounds {} {
	return "$Bounds(min,x) $Bounds(max,x) $Bounds(min,y) $Bounds(max,y) $Bounds(min,z) $Bounds(max,z)"
    }

    method ComputeBounds {} {
	if {$ComputeBoundsFlag} {

	    set ComputeBoundsFlag 0

	    foreach plane "x y z" { 
		# SGS where can we get these from, VTK class
		set Bounds(min,$plane) $::cvis::math::FLT_MAX
		set Bounds(max,$plane) [expr -$::cvis::math::FLT_MAX]
	    }

	    set collection $this.cvisActors
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {

		if {[$item GetVisibility]} {
		    scan [$item GetBounds] "%f %f %f %f %f %f" \
			    real(min,x) real(max,x) \
			    real(min,y) real(max,y) \
			    real(min,z) real(max,z)
		    
		    foreach plane "x y z" {
			if { $Bounds(min,$plane) > $real(min,$plane) } {
			    set Bounds(min,$plane) $real(min,$plane)
			}
			
			if { $Bounds(max,$plane) < $real(max,$plane) } {
			set Bounds(max,$plane) $real(max,$plane)
			}
		    }
		}
		
		set item [$collection GetNextItem]
	    }
	}
    }


    method Expose {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$InExpose == 1 } {
	    return 
	}
	
	set InExpose 1
	$RenderWindow SetDesiredUpdateRate $StillUpdateRate
	update
	RenderInteractive
	set InExpose 0
    }

    variable OldFocus ""

    method Enter {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]
	set OldFocus [focus]
	focus $widget
	UpdateRenderer $x $y
    }

    method Exit {} {
	set tracervar [::itcl::local cvisTrace #auto]
	# SGS seems to not work very well why to VTK folks want this
	# focus $OldFocus
    }

    method ExitApplication {} {
	::cvis::Exit
    }

    variable PickPointColor "1.0 1.0 1.0"

    method SetPickPointColor {r g b} {
	set PickPointColor "$r $g $b"
	if {[string length [info commands $this.pickerActor]]} {
	}

	if {[string length [info commands $this.boundaryActor]]} {
	    eval [$this.boundaryActor GetProperty] SetColor $PickPointColor
	}
    }

    variable PickPointIsVisible 1

    method SetPickPointVisibility {visible} {
	if {$visible != $PickPointIsVisible} {
	    set PickPointIsVisible $visible
	    if {[string compare $PickBounds "0 0 0 0 0 0"]} {
		$this.boundaryActor SetVisibility $PickPointIsVisible
	    } {
		$this.boundaryActor SetVisibility 0
	    }
	}
    }

    method GetPickPointVisibility {} {
	return $PickPointIsVisible
    }

    method SetPickPoint {x y z} {
	set PickPoint(x) $x
	set PickPoint(y) $y
	set PickPoint(z) $z

	# Force A redraw of the text widget
	#		Update
	
	# Picker in render window is in 3D space
	set camera [$Renderer GetActiveCamera]
	
		$PickPoint(x) $PickPoint(y) $PickPoint(z)

	# Set the size of the picker
	
	set distance [$camera GetDistance] 

	set view_angle [$camera GetViewAngle]
	set view_angle [expr $::cvis::math::DEGREESTORADIANS * $view_angle]
	set scale [$camera GetParallelScale]
	
	set width [expr $distance * tan($view_angle / 2.0) * 0.01]
	
	
	# Text message for pick location
		"($PickPoint(x), $PickPoint(y), $PickPoint(z))"

	RenderInteractive
    }

    variable BoundaryEdgeWidth 3
    variable PickBounds "0 0 0 0 0 0"
    method SetPickBounds {min_x max_x min_y max_y min_z max_z} {
	set PickBounds "$min_x $max_x $min_y $max_y $min_z $max_z"
	    cvisActor $this.boundaryActor
	    eval [$this.boundaryActor GetProperty] SetColor $PickPointColor
	    
	    [$this.boundaryActor GetProperty] SetAmbient  1
	    [$this.boundaryActor GetProperty] SetDiffuse  0
	    [$this.boundaryActor GetProperty] SetSpecular 0
	    [$this.boundaryActor GetProperty] SetLineWidth \
		    $BoundaryEdgeWidth
	    $Renderer AddActor [$this.boundaryActor GetActor]
	}
	if {[string compare $PickBounds "0 0 0 0 0 0"]} {
	    $this.boundaryActor SetVisibility $PickPointIsVisible
	} {
	    $this.boundaryActor SetVisibility 0
	}

	RenderInteractive 
    }

    method GetPickBounds {} {
	return $PickBounds
    }

    method PickMethod { } {
	} {
		    pick_position(x) pick_position(y) pick_position(z)
	    
	    SetPickPoint $pick_position(x) $pick_position(y) $pick_position(z)

	    foreach object $Pickers { 
		$object SetPickPoint \
			$pick_position(x) $pick_position(y) $pick_position(z)
	    }
	}
    }

    variable AutoResize 1
    method SetAutoResize {resize} {
	set AutoResize $resize
    }

    method GetAutoResize {} {
	return $AutoResize
    }

    method DecorationHeight {} {
	    
    }

    method DecorationWidth {} {
    }

    method ResizeSize {widget width height} {
	if {!$InSetSize} {
		if {$AutoResize} {
		    # Account for toolbar and button bar
		    set Size(x) [expr $Size(x) - [DecorationWidth]]
		    set Size(y) [expr $Size(y) - [DecorationHeight]]
		    if {[string length $Interface]} {
			$Interface SetSize $Size(x) $Size(y)
		    } {
			SetSize $Size(x) $Size(y)
		    }
		}
	    }
	}
    }

    variable Pickers ""
    method AddPicker {object} {
	lappend Pickers $object
    }

    constructor {name} {
	global cvis_library

	::cvis::RestartFile::Add $this

	set Size(x) 640 
	set Size(y) 480 

	set Name $name

	cvisView $this.views
	cvisCollection $this.cvisActors
	cvisCollection $this.vtkActors


	# Need to turn off the culler so we don't drop small items from
	# the render
	set cullers [$Renderer GetCullers]
	$cullers InitTraversal
	set culler [$cullers GetNextItem]
	$culler SetMinimumCoverage 1e-12

	# Objects for the rubber band zoom feature


	$RubberBandPolyData SetPoints      $RubberBandPoints
	$RubberBandPolyData SetLines       $RubberBandLines
	$RubberBandMapper   SetInput       $RubberBandPolyData
	$RubberBandMapper   SetLookupTable $RubberBandColors
	$RubberBandActor    SetMapper      $RubberBandMapper

	$RubberBandColors SetNumberOfTableValues 2
	$RubberBandColors SetNumberOfColors 2
	$RubberBandColors SetTableValue 0 1.0 0.0 0.0 1.0
	$RubberBandColors SetTableValue 1 1.0 1.0 1.0 1.0
	
	[$RubberBandPolyData GetPointData] SetScalars $RubberBandScalars
	
	$RubberBandMapper SetScalarRange 0 1

	$RubberBandPoints InsertPoint 0  0  0  0
	$RubberBandPoints InsertPoint 1  0 10  0
	$RubberBandPoints InsertPoint 2 10 10  0
	$RubberBandPoints InsertPoint 3 10  0  0
	
	$RubberBandLines  InsertNextCell 5
	$RubberBandLines  InsertCellPoint 0
	$RubberBandLines  InsertCellPoint 1
	$RubberBandLines  InsertCellPoint 2
	$RubberBandLines  InsertCellPoint 3
	$RubberBandLines  InsertCellPoint 0
	
	$RubberBandScalars InsertNextTuple1 0
	$RubberBandScalars InsertNextTuple1 1
	$RubberBandScalars InsertNextTuple1 0
	$RubberBandScalars InsertNextTuple1 1

	$RubberBandMapper ScalarVisibilityOn

	set PickedAssembly ""
	set PrePickedProperty ""

	
	
	
	
	
	
	
	iwidgets::menubar $frame.menuBar -menubuttons {
	    menubutton file -text "File" -menu {
		options -tearoff false

		command Open -label "Open" -command {$this Open} \
			-helpstr "Open saved state file"

		command SaveAs -label "SaveAs" -command {$this SaveAs} \
			-helpstr "Save state to a file"

		command exit -label "Exit" -command {$this ExitApplication} \
			-helpstr "Exit Vizamrai"
	    }
	    menubutton help -text "Help" -menu {
		options -tearoff false

		command help -label "Help Contents" \
			-command {$this HelpContents} \
			-helpstr "Help Contents"
		
		command about -label "About" -command {$this About} \
			-helpstr "About Vizamrai"
	    }
	}
	
	pack $frame.menuBar -side top -fill x

	#****************************************
	#****************************************
	# SGS Working here

	frame $frame.sidebar
	pack $frame.sidebar -side left -fill y

	frame $frame.sidebar.empty

	frame $frame.sidebar.tools

	set toolbar $frame.sidebar.tools.view
	iwidgets::toolbar $toolbar -orient vertical

	set widget [$toolbar add checkbutton center \
		-command "$this Reset" \
		-balloonstr "Center" \
		-image [image create photo -file \
		        [file join $cvis_library images resetoff.gif]] \
		-selectimage [image create photo -file \
  		        [file join $cvis_library images reseton.gif]] \
			]

	$widget configure \
		-borderwidth 0 \
		-activebackground [$widget cget -background]

	bind $widget <Button-1> "%W select"

	set widget [$toolbar add checkbutton pan \
  		-balloonstr "Pan" \
		-image [image create photo -file \
  		        [file join $cvis_library images panoff.gif]] \
		-selectimage [image create photo -file \
  		        [file join $cvis_library images panon.gif]] \
			]

	$widget configure \
		-borderwidth 0 \
		-activebackground [$widget cget -background]

	$widget configure -command "$widget deselect"

	bind $widget <B1-Motion> "$this Pan %x %y"
	bind $widget <Button-1> "$this StartPan %W %x %y"
	bind $widget <ButtonRelease-1> "$this EndPan %W %x %y"

	set widget [$toolbar add checkbutton rotate \
  		-balloonstr "Rotate" \
  		-image [image create photo -file \
  		        [file join $cvis_library images rotateoff.gif]] \
		-selectimage [image create photo -file \
  		        [file join $cvis_library images rotateon.gif]] \
			]

	$widget configure \
		-borderwidth 0 \
		-activebackground [$widget cget -background]


	$widget configure -command "$widget deselect"

	bind $widget <B1-Motion> "$this Rotate %x %y"
	bind $widget <ButtonPress-1> "$this StartRotate %W %x %y"
	bind $widget <ButtonRelease-1> "$this EndRotate %W %x %y"

	set widget [$toolbar add checkbutton zoom \
  		-balloonstr "Zoom" \
  		-image [image create photo -file \
  		        [file join $cvis_library images zoomoff.gif]] \
		-selectimage [image create photo -file \
  		        [file join $cvis_library images zoomon.gif]] \
			]

	$widget configure \
		-borderwidth 0 \
		-activebackground [$widget cget -background]


	$widget configure -command "$widget deselect"

	bind $widget <B1-Motion> "$this Zoom %x %y"
	bind $widget <ButtonPress-1> "$this StartZoom %W %x %y"
	bind $widget <ButtonRelease-1> "$this EndZoom %W %x %y"

	pack $toolbar -anchor n -side top -fill y
	# ******************************************



	set toolbar $frame.sidebar.tools.user
	iwidgets::toolbar $toolbar -orient vertical

	set widget [$toolbar add button birdseye \
		-command "$this BirdsEye" \
		-balloonstr "BirdsEye" \
		-image [image create photo -file \
		        [file join $cvis_library images bird.gif]]]
	$widget configure -bd 0

	frame $toolbar.memorydot
	iwidgets::toolbar $toolbar.memorydot.left -orient vertical
	foreach i {1 3 5} {
	    set button [$toolbar.memorydot.left add button b$i \
		    -command "$this GetView $i" \
		    -balloonstr "View $i" \
		-image [image create photo -file \
		        [file join $cvis_library images memorydot.gif]]]
	    $button configure -bd 0
	    bind $button <B3-ButtonRelease> "$this SetView $i"
	}

	iwidgets::toolbar $toolbar.memorydot.right -orient vertical
	foreach i {2 4 6} {
	    set button [$toolbar.memorydot.right add button b$i \
		    -command "$this GetView $i" \
		    -balloonstr "View $i" \
		-image [image create photo -file \
		        [file join $cvis_library images memorydot.gif]]]
	    $button configure -bd 0
	    bind $button <B3-ButtonRelease> "$this SetView $i"
	}

	pack $toolbar.memorydot.left -side left 
	pack $toolbar.memorydot.right -side right
	pack $toolbar.memorydot -side bottom -fill y

	pack $toolbar -side bottom -fill y

	# *************************************************

	set toolbar $frame.sidebar.tools.presetl
	iwidgets::toolbar $toolbar -orient vertical

	set widget [$toolbar add button top \
		-command "$this Top" \
		-balloonstr "Top" \
		-image [image create photo -file \
		        [file join $cvis_library images top.gif]]]
	$widget configure -bd 0

	set widget [$toolbar add button left \
		-command "$this Left" \
		-balloonstr "Left" \
		-image [image create photo -file \
		        [file join $cvis_library images left.gif]]]
	$widget configure -bd 0

	set widget [$toolbar add button front \
		-command "$this Front" \
		-balloonstr "Front" \
		-image [image create photo -file \
		        [file join $cvis_library images front.gif]]]
	$widget configure -bd 0

	pack $toolbar -anchor n -side left -fill y

	# **********************************************8

	set toolbar $frame.sidebar.tools.presetr
	iwidgets::toolbar $toolbar -orient vertical

	set widget [$toolbar add button bottom \
		-command "$this Bottom" \
		-balloonstr "Bottom" \
		-image [image create photo -file \
		        [file join $cvis_library images bottom.gif]]]
	$widget configure -bd 0

	set widget [$toolbar add button right \
		-command "$this Right" \
		-balloonstr "Right" \
		-image [image create photo -file \
		        [file join $cvis_library images right.gif]]]
	$widget configure -bd 0

	set widget [$toolbar add button back \
		-command "$this Back" \
		-balloonstr "Back" \
		-image [image create photo -file \
		        [file join $cvis_library images back.gif]]]
	$widget configure -bd 0

	pack $toolbar -anchor n -side left -fill y

	#****************************************
	#****************************************

	set widget [button $frame.sideToolbarControl  \
		-command "$this HideSideToolbar $frame.sidebar" \
		-image [image create photo -file \
		        [file join $cvis_library images verthandleoff.gif]]]

	$widget configure \
		-borderwidth 0 \
		-activebackground [$widget cget -background]

	pack $widget -side left -fill y

	pack $frame.sidebar.tools

	# SGS Working here
	#****************************************
	#****************************************

	frame $frame.bottombar
	pack $frame.bottombar -side bottom -fill x

	frame $frame.bottombar.empty
	frame $frame.bottombar.tools

	set toolbar $frame.bottombar.tools.view
	iwidgets::toolbar $toolbar

	$toolbar add button saveImage \
		-command "$this Write" \
		-balloonstr "Save image to file" \
		-image [image create photo -file \
		        [file join $cvis_library images save.gif]]

	$toolbar add button render \
		-command "$this RenderNow" \
		-balloonstr "Render" \
		-image 	[image create photo -file \
		[file join $cvis_library images render.gif]]

	set RenderPipelineOnImage [image create photo -file \
		[file join $cvis_library images pipeon.gif]]

	set RenderPipelineOffImage [image create photo -file \
		[file join $cvis_library images pipeoff.gif]]

	set RenderPipelineWidget [$toolbar add button pipe \
		-command "$this RenderPipeToggle" \
		-balloonstr "Render Pipeline" \
		-image $RenderPipelineOnImage]

	set RenderWireFrameImage [image create photo -file \
		[file join $cvis_library images wireframe.gif]]

	set RenderFullFrameImage [image create photo -file \
		[file join $cvis_library images fullframe.gif]]

	set RenderWireFrameWidget [$toolbar add button wireframe \
		-command "$this RenderWireFrameToggle" \
		-balloonstr "Wireframe/Full Render" \
		-image $RenderFullFrameImage]

	$toolbar add button command \
		-command "cvisTkCon" \
		-balloonstr "Command Line Interface" \
		-image [image create photo -file \
		        [file join $cvis_library images command.gif]]

	pack $toolbar -side bottom -fill x 

	set widget [button $frame.bottomToolbarControl  \
		-command "$this HideBottomToolbar $frame.bottombar" \
		-image [image create photo -file \
		        [file join $cvis_library images horizhandleoff.gif]]]

	$widget configure \
		-borderwidth 0 \
		-activebackground [$widget cget -background]

	pack $widget -side bottom -fill x 

	pack $frame.bottombar.tools -fill x 

	set RenderWidget $widget
	
	vtkTkRenderWidget $widget \
		-width $Size(x) -height $Size(y)
	pack $widget -expand true
	
	# Bindings
	bind $widget <Expose> "$this Expose"
	bind $widget <Enter> "$this Enter %W %x %y"
	bind $widget <Leave> "$this Exit"
	
	bind $widget <Any-ButtonPress> "$this StartMotion %W %x %y"
	bind $widget <Any-ButtonRelease> "$this EndMotion %W %x %y"
	bind $widget <B1-Motion> "$this Rotate %x %y"
	bind $widget <B2-Motion> "$this Pan %x %y"
	bind $widget <Shift-B1-Motion> "$this Pan %x %y"
	bind $widget <B3-Motion> "$this Zoom %x %y"
        bind $widget <Shift-B3-Motion> "$this RubberZoom %W %x %y"

	# This currently will only work on Windows only.
	bind $widget <MouseWheel> "$this WheelZoom %D"
	
	bind $widget <KeyPress-c> {cvisTkCon}
	bind $widget <KeyPress-r> "$this ResetInteractor %W %x %y"
	bind $widget <KeyPress-f> "$this ResetInteractorFullScreen %W %x %y"
	
	bind $widget <KeyPress-w> "$this Wireframe"
	bind $widget <KeyPress-s> "$this Surface"
	
	bind $widget <KeyPress-space> "$this PickActor %W %x %y"

	bind $widget <KeyPress-o> "$this Write"


	set RenderWindow [$widget GetRenderWindow]

	$RenderWindow AddRenderer $Renderer

	set WindowNameProp [cvisProperty $name.WindowName]
	$name.WindowName AddCallback "$this UpdateWindowName"

	# SGS this was added to create a default light but 
	# there should be some better way of doing this?
	# If this is not here then we get no light and nothing
	# ever shows up
	$Renderer Render
	UpdateRenderer 0 0 

	if { [file exists $::cvis::options::DEFAULT_CAMERA] } {
	    RestoreView	 $::cvis::options::DEFAULT_CAMERA
	}

	foreach plane "x y z" { 
	    set Bounds(min,$plane) 0.0
	    set Bounds(max,$plane) 1.0
	}

	foreach light {1 2 3 4} {
	    vtkLight $this.light$light
	    $this.light$light SetLightType 3
	    $this.light$light SetSwitch 0
	    $Renderer AddLight $this.light$light
	}

	
    }

    destructor {



  	$Renderer Delete
    }

    method GetImageWriter { } {
	return [namespace current]::$ImageWriter
    }
    
    method GetPicker { } {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method GetRenderWindowInteractor { } {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method GetRenderer { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Renderer
    }

    method GetWindowName {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$WindowNameProp GetValue]
    }

    method GetWindowNameProp {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$WindowNameProp GetReference]
    }

    method SetWindowName {windowname} {
	set tracervar [::itcl::local cvisTrace #auto]
	$WindowNameProp SetValue $windowname
    }
    
    method SetEquivWindowName { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$WindowNameProp SetEquiv $property
    }

    method UpdateWindowName { windowname } {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method GetRenderWindow { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $RenderWindow
    }

    method SetBackground { r g b } {
	set tracervar [::itcl::local cvisTrace #auto]
	$Renderer SetBackground $r $g $b
	Render
    }

    method GetBackground {r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Renderer GetBackground]
    }

    variable RenderWireFrameWidget 
    variable RenderWireFrame 0
    variable RenderWireFrameImage
    variable RenderFullFrameImage
    method RenderWireFrameToggle {} {
	if {$RenderWireFrame} {
	    set RenderWireFrame 0
	    $RenderWireFrameWidget configure -image $RenderFullFrameImage
	    Surface
	} {
	    set RenderWireFrame 1
	    $RenderWireFrameWidget configure -image $RenderWireFrameImage
	    Wireframe
	}

	$RenderWindow SetDesiredUpdateRate $StillUpdateRate
    }

    variable RenderPipelineWidget
    variable RenderPipeline 1
    variable RenderPipelineOnImage
    variable RenderPipelineOffImage
    method RenderPipeToggle {} {
	if {$RenderPipeline} {
	    set RenderPipeline 0
	    ::cvis::queue::SetPipeline 0
	    $RenderPipelineWidget configure -image $RenderPipelineOffImage

	} {
	    set RenderPipeline 1
	    ::cvis::queue::SetPipeline 1
	    $RenderPipelineWidget configure -image $RenderPipelineOnImage
	}
    }

    # Force a Render now, regardless of the current pipeline state
    method RenderNow {} {
	set tracervar [::itcl::local cvisTrace #auto]

	Render

	if {![::cvis::queue::GetPipeline]} {
	    ::cvis::queue::SetPipeline 1
	    ::cvis::queue::SetPipeline 0
	}
    }

    # Schedule a Render 
    method Render { } {
	set tracervar [::itcl::local cvisTrace #auto]
	
	# Schedule a render operation
	::cvis::queue::Add "$this Execute"
	
	# Excute the waiting operations, if any
	::cvis::queue::Execute
    }

    # Called by interactive changes (things which 
    # only effect the renderer, not other elements)
    method RenderInteractive {} {
	# If interactive just render now otherwise schedule
	# a render for later
	if {$RenderPipeline} {
	    $RenderWindow Render    
	} {
	    Render
	}
    }

    method Execute { } {
	set tracervar [::itcl::local cvisTrace #auto]

	$RenderWindow Render    
    }

    method Write {} {
	set tracervar [::itcl::local cvisTrace #auto]
	set filename [GetFileName]
	if {[string length [GetFileName]]} {
	    RenderNow
	}
    }

    method SetFileName { filename } {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method GetFileName { } {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    variable InSetSize 0

    method SetSize { size_x size_y } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Size(x) $size_x
	set Size(y) $size_y
	set InSetSize 1
	

	if { ($Size(x) > $current_max_width) || \
		($Size(y) > $current_max_height) } {
		    "[expr $Size(x) + [DecorationWidth]]x[expr $Size(y) + [DecorationHeight]]"

	} 

		-height $Size(y)
	update
	set InSetSize 0
    }

    method SaveState { filename } {
	set tracervar [::itcl::local cvisTrace #auto]
	SaveView [$Renderer GetActiveCamera] $filename
    }

    method RestoreState { filename} {
	set tracervar [::itcl::local cvisTrace #auto]
	RestoreView [$Renderer GetActiveCamera] $filename
	Render
    }

    method AddActor { actor } {
	set tracervar [::itcl::local cvisTrace #auto]
	if { ![string compare [$actor GetClassName] "cvisActor"] } {
	    $this.cvisActors AddItem $actor
	    $Renderer AddActor [$actor GetActor]
	} {
	    $this.vtkActors AddItem $actor
	    $Renderer AddActor $actor
	}
	set ComputeBoundsFlag 1
    }

    method RemoveActor { actor } {
	set tracervar [::itcl::local cvisTrace #auto]
	if { ![string compare [$actor GetClassName] "cvisActor"] } {
	    $Renderer RemoveActor [$actor GetActor]
	    $this.cvisActors RemoveItem $actor
	} {
	    $Renderer RemoveActor $actor
	    $this.vtkActors RemoveItem $actor
	}
	set ComputeBoundsFlag 1
    }

    method AddVolume { volume } {
	set tracervar [::itcl::local cvisTrace #auto]
	$Renderer AddVolume $volume
	set ComputeBoundsFlag 1
    }

    method RemoveVolume { volume } {
	set tracervar [::itcl::local cvisTrace #auto]
	$Renderer RemoveVolume $volume
	set ComputeBoundsFlag 1
    }

    method SetOffScreenRendering {value} {
	set tracervar [::itcl::local cvisTrace #auto]
	$RenderWindow SetOffScreenRendering $value
    }

    method GetOffScreenRendering {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$RenderWindow GetOffScreenRendering]
    }

    method SetUniformScaling {uniformscaling} {
	set tracervar [::itcl::local cvisTrace #auto]
	set UniformScaling $uniformscaling
	
	SetScaling
    }

    method SetScaling { } {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$UniformScaling } { 

	    # Reset scaling for bounds computation
	    # This does not seem to cause a speed problem but is
	    # there a better way to compute this
	    set collection $this.cvisActors
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		$item SetScale 1.0 1.0 1.0
		set item [$collection GetNextItem]
	    }
	    
	    # Scale display objects so they have the same size
	    # along each of their dimensions.  
	    foreach plane "x y z" { 
		# SGS where can we get these from, VTK class
		set bounds(min,$plane) $::cvis::math::FLT_MAX
		set bounds(max,$plane) [expr -$::cvis::math::FLT_MAX]
	    }

	    set collection $this.cvisActors
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		scan [$item GetBounds] "%f %f %f %f %f %f" \
			real(min,x) real(max,x) \
			real(min,y) real(max,y) \
			real(min,z) real(max,z)
		
		foreach plane "x y z" {
		    if { $bounds(min,$plane) > $real(min,$plane) } {
			set bounds(min,$plane) $real(min,$plane)
		    }
		    
		    if { $bounds(max,$plane) < $real(max,$plane) } {
			set bounds(max,$plane) $real(max,$plane)
		    }
		}
		
		set item [$collection GetNextItem]
	    }

	    # SGS instead of scaling to 1000.0 what about scaling to
	    # Max extent size to avoid size changing as much

	    set scaleto [expr -$::cvis::math::FLT_MAX]
	    foreach plane "x y z" { 
		if { $bounds(max,$plane) > $scaleto } {
		    set scaleto $bounds(max,$plane)
		}
	    }
		

	    set Scaling(x) [expr $scaleto / \
		    ( $bounds(max,x) - $bounds(min,x) )]
	    set Scaling(y) [expr $scaleto / \
		    ( $bounds(max,y) - $bounds(min,y) )]
	    set Scaling(z) [expr $scaleto / \
		    ( $bounds(max,y) - $bounds(min,y) )]
	} {
	    
	    # Do not do any scaling
	    set Scaling(x) 1.0
	    set Scaling(y) 1.0
	    set Scaling(z) 1.0
	}
       

	set collection "$this.cvisActors"
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    $item SetScale $Scaling(x) $Scaling(y) $Scaling(z)
	    set item [$collection GetNextItem]
	}
    }

    variable AutoSave 0

    method SetAutoSave { autosave } {
	set tracervar [::itcl::local cvisTrace #auto]
	set AutoSave $autosave
    }

    method GetAutoSave { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $AutoSave
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]
	# SGS If autosave turned on then save the image to the current
	# file
	if {$AutoSave} { 
	    ::cvis::queue::Add "$this Write"
	}
    }

    variable FileBrowser ".null"

    # Save the current state of the camera to a file.  
    # If not filename is passed in a file browser will be displayed
    method SaveView {{filename ""}} {
	set tracervar [::itcl::local cvisTrace #auto]

	set camera [[GetRenderer] GetActiveCamera]

	set types {
	    {{vtk Camera }  {.cam}        }
	    {{All Files }   *             }
	}
	
	if { $filename == "" } {
	    if {[string length [find objects $FileBrowser]]} {
		$FileBrowser SetSavePrompt
	    } {
		set FileBrowser [cvisFileBrowser #auto]
		$FileBrowser SetTitle "Select Camera File"
		$FileBrowser SetFileTypes $types
		$FileBrowser SetSavePrompt
	    }

	    set filename [$FileBrowser Prompt]
	}

	if { $filename != "" } {
	    set fp [open $filename w]
	    if { $fp == "" } {
		puts "Error opening camera view file"
	    } else {
		puts $fp "# vtk CameraFile Version 2.0"
		puts $fp "Clipping Range: [$camera GetClippingRange]"
		puts $fp "Distance: [$camera GetDistance]"
		puts $fp "Eye Angle: [$camera GetEyeAngle]"
		puts $fp "Focal Disk: [$camera GetFocalDisk]"
		puts $fp "Focal Point: [$camera GetFocalPoint]"
		puts $fp "Position: [$camera GetPosition]"
		puts $fp "Parallel Projection: [$camera GetParallelProjection]"
		puts $fp "Thickness: [$camera GetThickness]"
		puts $fp "View Angle: [$camera GetViewAngle]"
		puts $fp "View Up: [$camera GetViewUp]"
		close $fp
	    }
	}
    }

    ## Recover view from saved file
    ## Provide camera name
    method RestoreView {{filename ""}} {
	set tracervar [::itcl::local cvisTrace #auto]

	set camera [[GetRenderer] GetActiveCamera]

	set types {
	    {{vtk Camera }  {.cam}        }
	    {{All Files }   *             }
	}

	if { $filename == "" } {
	    if {[string length [find objects $FileBrowser]]} {
		$FileBrowser SetOpenPrompt
	    } {
		set FileBrowser [cvisFileBrowser #auto]
		$FileBrowser SetTitle "Select Camera File"
		$FileBrowser SetFileTypes $types
		$FileBrowser SetOpenPrompt
	    }

	    set filename [$FileBrowser Prompt]
	}
	
	if { $filename != "" } {
	    set fp [open $filename r]
	    if { $fp == "" } {
		puts "Error opening camera view file"
	    } else {
		set line [read $fp]
		
		if { [string first "# vtk CameraFile Version 2.0" $line] == -1} {
		    puts "Can't read this file, can't find cookie"
		} else {
		    set loc [string first "Clipping Range:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 15] end]
			scan $substr "%f %f" min max
			$camera SetClippingRange $min $max
		    }
		    
		    set loc [string first "Distance:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 9] end]
			scan $substr "%f" distance
			$camera SetDistance $distance
		    }
		    
		    set loc [string first "Eye Angle:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 10] end]
			scan $substr "%f" angle
			$camera SetEyeAngle $angle
		    }
		    
		    set loc [string first "Focal Disk:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 11] end]
			scan $substr "%f" fd
			$camera SetFocalDisk $fd
		    }
		    
		    set loc [string first "Focal Point:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 12] end]
			scan $substr "%f %f %f" fp1 fp2 fp3
			$camera SetFocalPoint $fp1 $fp2 $fp3
		    }
		    
		    set loc [string first "Position:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 9] end]
			scan $substr "%f %f %f" p1 p2 p3
			$camera SetPosition $p1 $p2 $p3
		    }
		    
		    set loc [string first "Parallel Projection:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 20] end]
			scan $substr "%d" pp
			$camera SetParallelProjection $pp
		    }
		    
		    set loc [string first "Thickness:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 10] end]
			scan $substr "%f" t
			$camera SetThickness $t
		    }
		    
		    set loc [string first "View Angle:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 11] end]
			scan $substr "%f" va
			$camera SetViewAngle $va
		    }
		    
		    set loc [string first "View Up:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 8] end]
			scan $substr "%f %f %f" vup1 vup2 vup3
			$camera SetViewUp $vup1 $vup2 $vup3
		    }
		}

		# Need to force a render otherwise will only
		# happen when user causes a window update
		[GetRenderWindow] Render
	    }
	}
    }

    variable FileName ""
    variable FileTypes {{{Vizamrai State File} {.vstate}}}

    method Open {} {
	set ImageFileBrowser [cvisFileBrowser #auto]
	$ImageFileBrowser SetTitle "Select Vizamrai Save State File"
	$ImageFileBrowser SetFileTypes $FileTypes
	$ImageFileBrowser SetOpenPrompt
	set FileName [$ImageFileBrowser Prompt]
	::cvis::RestartFile::Open $FileName
    }

    method SaveAs {} {
	set ImageFileBrowser [cvisFileBrowser #auto]
	$ImageFileBrowser SetTitle "Select Vizamrai Save State File"
	$ImageFileBrowser SetFileTypes $FileTypes
	$ImageFileBrowser SetSavePrompt
	set FileName [$ImageFileBrowser Prompt]
	::cvis::RestartFile::Save $FileName
    }

    method HelpContents { } {
	set tracervar [::itcl::local cvisTrace #auto]
	global cvis_library 
	if { [winfo exists .help] } {
	    .help activate 
	} {
	    iwidgets::hyperhelp .help \
		    -title "Vizamrai Help" \
		    -modality none \
		    -topics {Contents} \
		    -helpdir [file join $cvis_library help]
	    .help showtopic Contents
	    .help activate
	}
    }

    method About { } {
	set tracervar [::itcl::local cvisTrace #auto]
	::cvis::gui::About
    }

    # SGS modified from DoRubberZoom
    method DoFullZoom { widget depth factor} {
	set tracervar [::itcl::local cvisTrace #auto]

	# Return if there is no renderer, or the rubber band is less
	# that 5 pixels in either direction
	if { ! $RendererFound } { return }
	if { [expr $StartRubberZoomX - $EndRubberZoomX] < 5 && \
		[expr $StartRubberZoomX - $EndRubberZoomX] > -5 } { return }
	if { [expr $StartRubberZoomY - $EndRubberZoomY] < 5 && \
		[expr $StartRubberZoomY - $EndRubberZoomY] > -5 } { return }
	
	# We'll need the window height later
	set WindowY [expr double([lindex [$widget configure -height] 4])]
	
	# What is the center of the rubber band box in pixels?
	set centerX [expr ($StartRubberZoomX + $EndRubberZoomX)/2.0]
	set centerY [expr ($StartRubberZoomY + $EndRubberZoomY)/2.0]
	
	# Convert the focal point to a display coordinate in order to get the
	# depth of the focal point in display units
	set FPoint [$CurrentCamera GetFocalPoint]
        set FPoint0 [expr double([lindex $FPoint 0])]
        set FPoint1 [expr double([lindex $FPoint 1])]
        set FPoint2 [expr double([lindex $FPoint 2])]
	$Renderer SetWorldPoint $FPoint0 $FPoint1 $FPoint2 1.0
	$Renderer WorldToDisplay
	set DPoint [$Renderer GetDisplayPoint]
	set focalDepth [expr double([lindex $DPoint 2])]
	
	# Convert the position of the camera to a display coordinate in order
	# to get the depth of the camera in display coordinates. Note this is
	# a negative number (behind the near clipping plane of 0) but it works
	# ok anyway
	set PPoint [$CurrentCamera GetPosition]
        set PPoint0 [expr double([lindex $PPoint 0])]
        set PPoint1 [expr double([lindex $PPoint 1])]
        set PPoint2 [expr double([lindex $PPoint 2])]
	$Renderer SetWorldPoint $PPoint0 $PPoint1 $PPoint2 1.0
	$Renderer WorldToDisplay
	set DPoint [$Renderer GetDisplayPoint]
	set positionDepth [expr double([lindex $DPoint 2])]
	
	# Find out the world position of where our new focal point should
	# be - it will be at the center of the box, back at the same focal depth
	# Don't actually set it now - we need to do all our computations before
	# we modify the camera
	$Renderer SetDisplayPoint $centerX $centerY $focalDepth
	$Renderer DisplayToWorld
	set newFocalPoint [$Renderer GetWorldPoint]
	set newFocalPoint0 [expr double([lindex $newFocalPoint 0])]
	set newFocalPoint1 [expr double([lindex $newFocalPoint 1])]
	set newFocalPoint2 [expr double([lindex $newFocalPoint 2])]
	set newFocalPoint3 [expr double([lindex $newFocalPoint 3])]
	if { $newFocalPoint3 != 0.0 } {
	    set newFocalPoint0 [expr $newFocalPoint0 / $newFocalPoint3]
	    set newFocalPoint1 [expr $newFocalPoint1 / $newFocalPoint3]
	    set newFocalPoint2 [expr $newFocalPoint2 / $newFocalPoint3]
	}
	
	# Find out where the new camera position will be - at the center of
	# the rubber band box at the position depth. Don't set it yet...
	$Renderer SetDisplayPoint $centerX $centerY $positionDepth
	$Renderer DisplayToWorld
	set newPosition [$Renderer GetWorldPoint]
	set newPosition0 [expr double([lindex $newPosition 0])]
	set newPosition1 [expr double([lindex $newPosition 1])]
	set newPosition2 [expr double([lindex $newPosition 2])]
	set newPosition3 [expr double([lindex $newPosition 3])]
	if { $newPosition3 != 0.0 } {
	    set newPosition0 [expr $newPosition0 / $newPosition3]
	    set newPosition1 [expr $newPosition1 / $newPosition3]
	    set newPosition2 [expr $newPosition2 / $newPosition3]
	}
	
	# We figured out how to position the camera to be centered, now we
	# need to "zoom". In parallel, this is simple since we only need to
	# change our parallel scale to encompass the entire y range of the
	# rubber band box. In perspective, we assume the box is drawn on the
	# near plane - this means that it is not possible that someone can
	# draw a rubber band box around a nearby object and dolly past it. It 
	# also means that you won't get very close to distance objects - but that
	# seems better than getting lost.
	if {[$CurrentCamera GetParallelProjection]} {
	    # the new scale is just based on the y size of the rubber band box
	    # compared to the y size of the window
	    set ydiff [expr ($StartRubberZoomX - $EndRubberZoomX)]
	    if { $ydiff < 0.0 } { set ydiff [expr $ydiff * -1.0] }
	    set newScale [$CurrentCamera GetParallelScale]
	    set newScale [expr $newScale * $ydiff / $WindowY]
	    
	    # now we can actually modify the camera
	    $CurrentCamera SetFocalPoint $newFocalPoint0 $newFocalPoint1 $newFocalPoint2
	    $CurrentCamera SetPosition $newPosition0 $newPosition1 $newPosition2
	    $CurrentCamera SetParallelScale $newScale
	    
	} else {
	    # find out the center of the rubber band box on the near plane
	    $Renderer SetDisplayPoint $centerX $centerY $depth
	    $Renderer DisplayToWorld
	    set nearFocalPoint [$Renderer GetWorldPoint]
	    set nearFocalPoint0 [expr double([lindex $nearFocalPoint 0])]
	    set nearFocalPoint1 [expr double([lindex $nearFocalPoint 1])]
	    set nearFocalPoint2 [expr double([lindex $nearFocalPoint 2])]
	    set nearFocalPoint3 [expr double([lindex $nearFocalPoint 3])]
	    if { $nearFocalPoint3 != 0.0 } {
		set nearFocalPoint0 [expr $nearFocalPoint0 / $nearFocalPoint3]
		set nearFocalPoint1 [expr $nearFocalPoint1 / $nearFocalPoint3]
		set nearFocalPoint2 [expr $nearFocalPoint2 / $nearFocalPoint3]
	    }
	    
	    # find the world coordinates of the point centered on the rubber band box
	    # in x, on the border in y, and at the near plane depth.
	    $Renderer SetDisplayPoint $centerX $StartRubberZoomY  $depth
	    $Renderer DisplayToWorld
	    set focalEdge [$Renderer GetWorldPoint]
	    set focalEdge0 [expr double([lindex $focalEdge 0])]
	    set focalEdge1 [expr double([lindex $focalEdge 1])]
	    set focalEdge2 [expr double([lindex $focalEdge 2])]
	    set focalEdge3 [expr double([lindex $focalEdge 3])]
	    if { $focalEdge3 != 0.0 } {
		set focalEdge0 [expr $focalEdge0 / $focalEdge3]
		set focalEdge1 [expr $focalEdge1 / $focalEdge3]
		set focalEdge2 [expr $focalEdge2 / $focalEdge3]
	    }
	    
	    # how far is this "rubberband edge point" from the focal point?
	    set ydist [expr \
		    sqrt( \
		    ($nearFocalPoint0 - $focalEdge0)*($nearFocalPoint0 - $focalEdge0) + \
		    ($nearFocalPoint1 - $focalEdge1)*($nearFocalPoint1 - $focalEdge1) + \
		    ($nearFocalPoint2 - $focalEdge2)*($nearFocalPoint2 - $focalEdge2) )]
	    
	    # We need to know how far back we must be so that when we view the scene
	    # with the current view angle, we see all of the y range of the rubber
	    # band box. Use a simple tangent equation - opposite / adjacent = tan theta
	    # where opposite is half the y height of the rubber band box on the near
	    # plane, adjacent is the distance we are solving for, and theta is half
	    # the viewing angle. This distance that we solve for is the new distance
	    # to the near plane - to find the new distance to the focal plane we
	    # must take the old distance to the focal plane, subtract the near plane
	    # distance, and add in the distance we solved for.
	    set angle [expr 0.5 * (3.141592 / 180.0) * [$CurrentCamera GetViewAngle]]
	    set d [expr $ydist/tan($angle)]
	    set range [$CurrentCamera GetClippingRange]
	    set nearplane [expr double([lindex $range 0])]

	    # now we can actually modify the camera
	    $CurrentCamera SetFocalPoint $newFocalPoint0 $newFocalPoint1 $newFocalPoint2
	    $CurrentCamera SetPosition $newPosition0 $newPosition1 $newPosition2
	    $CurrentCamera Dolly $factor
	    $Renderer ResetCameraClippingRange
	}    
    }

    method DoRubberZoom { widget } {
	set tracervar [::itcl::local cvisTrace #auto]

	# Return if there is no renderer, or the rubber band is less
	# that 5 pixels in either direction
	if { ! $RendererFound } { return }
	if { [expr $StartRubberZoomX - $EndRubberZoomX] < 5 && \
		[expr $StartRubberZoomX - $EndRubberZoomX] > -5 } { return }
	if { [expr $StartRubberZoomY - $EndRubberZoomY] < 5 && \
		[expr $StartRubberZoomY - $EndRubberZoomY] > -5 } { return }
	
	# We'll need the window height later
	set WindowY [expr double([lindex [$widget configure -height] 4])]
	
	# What is the center of the rubber band box in pixels?
	set centerX [expr ($StartRubberZoomX + $EndRubberZoomX)/2.0]
	set centerY [expr ($StartRubberZoomY + $EndRubberZoomY)/2.0]
	
	# Convert the focal point to a display coordinate in order to get the
	# depth of the focal point in display units
	set FPoint [$CurrentCamera GetFocalPoint]
        set FPoint0 [expr double([lindex $FPoint 0])]
        set FPoint1 [expr double([lindex $FPoint 1])]
        set FPoint2 [expr double([lindex $FPoint 2])]
	$Renderer SetWorldPoint $FPoint0 $FPoint1 $FPoint2 1.0
	$Renderer WorldToDisplay
	set DPoint [$Renderer GetDisplayPoint]
	set focalDepth [expr double([lindex $DPoint 2])]
	
	# Convert the position of the camera to a display coordinate in order
	# to get the depth of the camera in display coordinates. Note this is
	# a negative number (behind the near clipping plane of 0) but it works
	# ok anyway
	set PPoint [$CurrentCamera GetPosition]
        set PPoint0 [expr double([lindex $PPoint 0])]
        set PPoint1 [expr double([lindex $PPoint 1])]
        set PPoint2 [expr double([lindex $PPoint 2])]
	$Renderer SetWorldPoint $PPoint0 $PPoint1 $PPoint2 1.0
	$Renderer WorldToDisplay
	set DPoint [$Renderer GetDisplayPoint]
	set positionDepth [expr double([lindex $DPoint 2])]
	
	# Find out the world position of where our new focal point should
	# be - it will be at the center of the box, back at the same focal depth
	# Don't actually set it now - we need to do all our computations before
	# we modify the camera
	$Renderer SetDisplayPoint $centerX $centerY $focalDepth
	$Renderer DisplayToWorld
	set newFocalPoint [$Renderer GetWorldPoint]
	set newFocalPoint0 [expr double([lindex $newFocalPoint 0])]
	set newFocalPoint1 [expr double([lindex $newFocalPoint 1])]
	set newFocalPoint2 [expr double([lindex $newFocalPoint 2])]
	set newFocalPoint3 [expr double([lindex $newFocalPoint 3])]
	if { $newFocalPoint3 != 0.0 } {
	    set newFocalPoint0 [expr $newFocalPoint0 / $newFocalPoint3]
	    set newFocalPoint1 [expr $newFocalPoint1 / $newFocalPoint3]
	    set newFocalPoint2 [expr $newFocalPoint2 / $newFocalPoint3]
	}
	
	# Find out where the new camera position will be - at the center of
	# the rubber band box at the position depth. Don't set it yet...
	$Renderer SetDisplayPoint $centerX $centerY $positionDepth
	$Renderer DisplayToWorld
	set newPosition [$Renderer GetWorldPoint]
	set newPosition0 [expr double([lindex $newPosition 0])]
	set newPosition1 [expr double([lindex $newPosition 1])]
	set newPosition2 [expr double([lindex $newPosition 2])]
	set newPosition3 [expr double([lindex $newPosition 3])]
	if { $newPosition3 != 0.0 } {
	    set newPosition0 [expr $newPosition0 / $newPosition3]
	    set newPosition1 [expr $newPosition1 / $newPosition3]
	    set newPosition2 [expr $newPosition2 / $newPosition3]
	}
	
	# We figured out how to position the camera to be centered, now we
	# need to "zoom". In parallel, this is simple since we only need to
	# change our parallel scale to encompass the entire y range of the
	# rubber band box. In perspective, we assume the box is drawn on the
	# near plane - this means that it is not possible that someone can
	# draw a rubber band box around a nearby object and dolly past it. It 
	# also means that you won't get very close to distance objects - but that
	# seems better than getting lost.
	if {[$CurrentCamera GetParallelProjection]} {
	    # the new scale is just based on the y size of the rubber band box
	    # compared to the y size of the window
	    set ydiff [expr ($StartRubberZoomX - $EndRubberZoomX)]
	    if { $ydiff < 0.0 } { set ydiff [expr $ydiff * -1.0] }
	    set newScale [$CurrentCamera GetParallelScale]
	    set newScale [expr $newScale * $ydiff / $WindowY]
	    
	    # now we can actually modify the camera
	    $CurrentCamera SetFocalPoint $newFocalPoint0 $newFocalPoint1 $newFocalPoint2
	    $CurrentCamera SetPosition $newPosition0 $newPosition1 $newPosition2
	    $CurrentCamera SetParallelScale $newScale
	    
	} else {
	    # find out the center of the rubber band box on the near plane
	    $Renderer SetDisplayPoint $centerX $centerY 0.0
	    $Renderer DisplayToWorld
	    set nearFocalPoint [$Renderer GetWorldPoint]
	    set nearFocalPoint0 [expr double([lindex $nearFocalPoint 0])]
	    set nearFocalPoint1 [expr double([lindex $nearFocalPoint 1])]
	    set nearFocalPoint2 [expr double([lindex $nearFocalPoint 2])]
	    set nearFocalPoint3 [expr double([lindex $nearFocalPoint 3])]
	    if { $nearFocalPoint3 != 0.0 } {
		set nearFocalPoint0 [expr $nearFocalPoint0 / $nearFocalPoint3]
		set nearFocalPoint1 [expr $nearFocalPoint1 / $nearFocalPoint3]
		set nearFocalPoint2 [expr $nearFocalPoint2 / $nearFocalPoint3]
	    }
	    
	    # find the world coordinates of the point centered on the rubber band box
	    # in x, on the border in y, and at the near plane depth.
	    $Renderer SetDisplayPoint $centerX $StartRubberZoomY  0.0
	    $Renderer DisplayToWorld
	    set focalEdge [$Renderer GetWorldPoint]
	    set focalEdge0 [expr double([lindex $focalEdge 0])]
	    set focalEdge1 [expr double([lindex $focalEdge 1])]
	    set focalEdge2 [expr double([lindex $focalEdge 2])]
	    set focalEdge3 [expr double([lindex $focalEdge 3])]
	    if { $focalEdge3 != 0.0 } {
		set focalEdge0 [expr $focalEdge0 / $focalEdge3]
		set focalEdge1 [expr $focalEdge1 / $focalEdge3]
		set focalEdge2 [expr $focalEdge2 / $focalEdge3]
	    }
	    
	    # how far is this "rubberband edge point" from the focal point?
	    set ydist [expr \
		    sqrt( \
		    ($nearFocalPoint0 - $focalEdge0)*($nearFocalPoint0 - $focalEdge0) + \
		    ($nearFocalPoint1 - $focalEdge1)*($nearFocalPoint1 - $focalEdge1) + \
		    ($nearFocalPoint2 - $focalEdge2)*($nearFocalPoint2 - $focalEdge2) )]
	    
	    # We need to know how far back we must be so that when we view the scene
	    # with the current view angle, we see all of the y range of the rubber
	    # band box. Use a simple tangent equation - opposite / adjacent = tan theta
	    # where opposite is half the y height of the rubber band box on the near
	    # plane, adjacent is the distance we are solving for, and theta is half
	    # the viewing angle. This distance that we solve for is the new distance
	    # to the near plane - to find the new distance to the focal plane we
	    # must take the old distance to the focal plane, subtract the near plane
	    # distance, and add in the distance we solved for.
	    set angle [expr 0.5 * (3.141592 / 180.0) * [$CurrentCamera GetViewAngle]]
	    set d [expr $ydist/tan($angle)]
	    set range [$CurrentCamera GetClippingRange]
	    set nearplane [expr double([lindex $range 0])]
	    set factor [expr [$CurrentCamera GetDistance] / \
		    ([$CurrentCamera GetDistance] - $nearplane + $d)]
	    
	    # now we can actually modify the camera
	    $CurrentCamera SetFocalPoint $newFocalPoint0 $newFocalPoint1 $newFocalPoint2
	    $CurrentCamera SetPosition $newPosition0 $newPosition1 $newPosition2
	    $CurrentCamera Dolly $factor
	    $Renderer ResetCameraClippingRange
	}    
    }

    method RubberZoom {widget x y} {
	set tracervar [::itcl::local cvisTrace #auto]

	if { ! $RendererFound } { return }
	
	set WindowY [lindex [$widget configure -height] 4]
	
	if { ! $RubberZoomPerformed } {
	    $Renderer AddProp $RubberBandActor
	    
	    set StartRubberZoomX $x
	    set StartRubberZoomY [expr $WindowY - $y - 1]
	    
	    set RubberZoomPerformed 1
	}
	
	set EndRubberZoomX $x
	set EndRubberZoomY [expr $WindowY - $y - 1]
	
	$RubberBandPoints SetPoint 0 $StartRubberZoomX $StartRubberZoomY  0
	$RubberBandPoints SetPoint 1 $StartRubberZoomX $EndRubberZoomY    0
	$RubberBandPoints SetPoint 2 $EndRubberZoomX   $EndRubberZoomY    0
	$RubberBandPoints SetPoint 3 $EndRubberZoomX   $StartRubberZoomY  0
	
	Render
    }

    method GetLight {light} {
	if {[string compare "camera" $light]} {
	    return $this.light$light
	} {
	    # Camera light is the default camera, first in list
	    set lights [$Renderer GetLights]
	    $lights InitTraversal
	    return [$lights GetNextItem]
	}
    }
}







