##
## File:        cvis2DRenderWindow.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: 2D Render Window
##


class cvis2DRenderWindow {
    variable Name
    variable Create 1
    
    variable UniformScaling 0

    variable WindowNameProp

    variable Renderer ""
    variable RenderWindow ""

    variable Size

    variable Bounds 
    variable Delta 
    variable Scaling
    
    variable GlobalScaling 1.0

    # Interactors
    variable InExpose 0

    variable RendererFound 0

    variable LastX 0
    variable LastY 0

    variable WindowCenterX 0 
    variable WindowCenterY 0

    method UpdateRenderer {widget x y} {
	
	# Get the renderer window dimensions
	set WindowX [lindex [$widget configure -width] 4]
	set WindowY [lindex [$widget configure -height] 4]

	# Find which renderer event has occurred in
	set renderers [$RenderWindow GetRenderers]
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
		set WindowCenterX [expr double($WindowX)*(($vpxmax - $vpxmin)/2.0\
			+ $vpxmin)]
		set WindowCenterY [expr double($WindowY)*(($vpymax - $vpymin)/2.0\
			+ $vpymin)]
		break
	    }
	}

	set CurrentCamera [$Renderer GetActiveCamera]

	# SGS First light is "headlight"
	# This should probably be an explict variable in class 
	# This is from VTK stuff
	set lights [$Renderer GetLights]
	$lights InitTraversal; set CurrentLight [$lights GetNextItem]

	set LastX $x
	set LastY $y
    }

    variable CurrentCamera ""
    variable CurrentLight ""

    method StartMotion {widget x y} {

	UpdateRenderer $widget $x $y

	if { $RendererFound } { 
	    $RenderWindow SetDesiredUpdateRate 1.0
	}
    }

    method EndMotion {widget x y} {
	if { $RendererFound } {
	    $RenderWindow SetDesiredUpdateRate 0.01
	    Render
	}
    }

    method Wireframe {widget} {

	set actors [$Renderer GetActors]

	$actors InitTraversal
	set actor [$actors GetNextItem]
	while { $actor != "" } {
	    [$actor GetProperty] SetRepresentationToWireframe
	    set actor [$actors GetNextItem]
	}
	Render
    }

    method Surface {widget} {
	
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

	Render
    }

    method Reset {widget x y} {

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


    method Rotate {widget x y} {

	if { $RendererFound } {

	    $CurrentCamera Azimuth [expr ($LastX - $x)]
	    $CurrentCamera Elevation [expr ($y - $LastY)]
	    $CurrentCamera OrthogonalizeViewUp

	    set LastX $x
	    set LastY $y

	    Render 
	}
    }

    method Pan {widget x y} {

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

	# SGS Added this
	set YCenter [expr ($FPoint1 - $RPoint1)/2.0 + $FPoint1]
	set XCenter [expr ($FPoint0 - $RPoint0)/2.0 + $FPoint0] 
	
	$CurrentCamera SetPosition \
		[expr ($FPoint0 - $RPoint0)/2.0 + $PPoint0] \
		[expr ($FPoint1 - $RPoint1)/2.0 + $PPoint1] \
		[expr ($FPoint2 - $RPoint2)/2.0 + $PPoint2]
	
	set LastX $x
	set LastY $y
	
	Render 
    }

    method Zoom {widget x y} {
	
	if { ! $RendererFound } {
	    return 
	}

	set zoomFactor [expr pow(1.02,(0.5*($y - $LastY)))]
	
	if {[$CurrentCamera GetParallelProjection]} {
	    set parallelScale [expr [$CurrentCamera GetParallelScale] * $zoomFactor];
	    $CurrentCamera SetParallelScale $parallelScale;
	    # SGS Added
	    set Zoom [expr (1.0 / $parallelScale) * $BaseZoom]
	} else {
	    set clippingRange [$CurrentCamera GetClippingRange]
	    set minRange [lindex $clippingRange 0]
	    set maxRange [lindex $clippingRange 1]
	    $CurrentCamera SetClippingRange [expr $minRange / $zoomFactor] \
		    [expr $maxRange / $zoomFactor]
	    $CurrentCamera Dolly $zoomFactor
	}
	
	set LastX $x
	set LastY $y
	
	Render 
    }

    method Expose {} {
	
	if {$InExpose == 1 } {
	    return 
	}
	
	set InExpose 1
	update
	Render
	set InExpose 0
    }

    variable OldFocus ""

    method Enter {widget x y} {
	set OldFocus [focus]
	focus $widget
	UpdateRenderer $widget $x $y
    }

    method Exit {} {
	# SGS seems to not work very well why to VTK folks want this
	# focus $OldFocus
    }

    variable Plane "z"
    method SetPlane {plane} {
	set Plane $plane

	set camera [$Renderer GetActiveCamera]
	
	switch $Plane {
	    "x" {
		$camera SetViewUp   0 1 0
	    }
	    
	    "y" {
		$camera SetViewUp   0 0 1
	    }
	    "z" {
		$camera SetViewUp   0 1 0
	    }
	}
	$camera ComputeViewPlaneNormal
	Render
    }

    method SetCenter { x y } {

	set XCenter $x
	set YCenter $y

	SetYCenter
	SetXCenter
    }

    variable YCenter 0.0
    method SetYCenter { } {
	set CurrentCamera [$Renderer GetActiveCamera]
	
	set camera [$Renderer GetActiveCamera]
	scan [$camera GetFocalPoint] "%f %f %f" \
		focal_point(x) focal_point(y) focal_point(z)
	scan [$camera GetPosition] "%f %f %f" \
		position(x) position(y) position(z)

	# Based on the plane we are viewing set the camera
	# in the databrowser window

	switch $Plane { 
	    "x" {
		$camera SetPosition \
			1.0 $YCenter $position(z)
		$camera SetFocalPoint \
			0.0 $YCenter $focal_point(z)
	    }
	    
	    "y" {
		$camera SetPosition \
			$position(x) -1.0 $YCenter 
		$camera SetFocalPoint \
			$focal_point(x) 0.0 $YCenter 
	    }
	    
	    "z" {
		$camera SetPosition \
			$position(x) $YCenter 1.0
		$camera SetFocalPoint \
			$focal_point(x) $YCenter 0.0
	    }
	}
	$camera ComputeViewPlaneNormal

	Render 
    }

    variable XCenter 0.0
    method SetXCenter { } {
	set CurrentCamera [$Renderer GetActiveCamera]
	
	set camera [$Renderer GetActiveCamera]
	scan [$camera GetFocalPoint] "%f %f %f" \
		focal_point(x) focal_point(y) focal_point(z)
	scan [$camera GetPosition] "%f %f %f" \
		position(x) position(y) position(z)

	# Based on the plane we are viewing set the camera
	# in the databrowser window
	switch $Plane { 
	    "x" {
		$camera SetPosition \
			1.0 $position(y) $XCenter 
		$camera SetFocalPoint \
			0.0 $focal_point(y) $XCenter 
	    }
	    
	    "y" {
		$camera SetPosition \
			$XCenter -1.0 $position(z)
		$camera SetFocalPoint \
			$XCenter 0.0 $focal_point(z)
	    }
	    
	    "z" {
		$camera SetPosition \
			$XCenter $position(y) 1.0
		$camera SetFocalPoint \
			$XCenter $focal_point(y) 0.0
	    }
	}
	$camera ComputeViewPlaneNormal

	Render 
    }

    variable Zoom 1.0
    variable BaseZoom 1.0
    method SetZoomFactor { } {
	set camera [$Renderer GetActiveCamera]
	$camera SetParallelScale [expr (1.0 / $Zoom) * $BaseZoom]
	Render 
    }

    constructor {name} {

	# Set Default values
	set Size(x) 640
	set Size(y) 480

	foreach plane "x y z" { 
	    set Delta($plane)      1.0
	    set Bounds(min,$plane) 0.0
	    set Bounds(max,$plane) 1.0
	}

	set Name $name



	set PickedAssembly ""
	set PrePickedProperty ""
	
	
	iwidgets::menubar $frame.menuBar -menubuttons {
	    menubutton file -text "File" -menu {
		options -tearoff false
		
		command exit -label "Exit" -command {exit} \
			-helpstr "Exit Vizamrai"
	    }
	    menubutton help -text "Help" -menu {
		options -tearoff false
		
		command about -label "About" -command {$this About} \
			-helpstr "About Vizamrai"
	    }
	}
	
	scale $frame.yscroll \
		-orient vertical \
		-from $Bounds(min,y) -to $Bounds(max,y) \
		-resolution $Delta(y) \
		-digits 9 \
		-variable [scope YCenter] \
		-showvalue true

	bind $frame.yscroll <ButtonRelease> \
		"$this SetYCenter"

	label $frame.ylabel -text "$YLabel"

	scale $frame.xscroll \
		-orient horizontal \
		-from $Bounds(min,x) -to $Bounds(max,x) \
		-resolution $Delta(x) \
		-digits 9 \
		-variable [scope XCenter] \
		-showvalue true

	bind $frame.xscroll <ButtonRelease> \
		"$this SetXCenter"

	label $frame.xlabel -text "$XLabel"

	scale $frame.zoom \
		-orient vertical \
		-from 1 -to 10 \
		-resolution 1 \
		-digits 9 \
		-variable [scope Zoom] \
		-showvalue true

	bind $frame.zoom <ButtonRelease> \
		"$this SetZoomFactor"

	label $frame.zlabel -text "$ZLabel"
    
	set RenderWidget $widget
	
	vtkTkRenderWidget $widget \
		-width $Size(x) -height $Size(y)

	grid $frame.menuBar -row 0 -column 0 -columnspan 3 -sticky ew
	grid $frame.ylabel -row 1 -column 0 -sticky s
	grid $frame.zlabel -row 1 -column 2 -sticky s
	grid $frame.yscroll -row 1 -column 0 -sticky ns
	grid $widget -row 1 -rowspan 1 -column 1 -sticky news
	grid $frame.zoom -row 1 -column 2 -sticky ns
	grid $frame.xscroll -row 2 -column 1 -sticky ew
	grid $frame.xlabel -row 2 -column 0 -sticky e
	
	# Bindings
	bind $widget <Expose> "$this Expose"
	bind $widget <Enter> "$this Enter %W %x %y"
	bind $widget <Leave> "$this Exit"
	
	bind $widget <Any-ButtonPress> "$this StartMotion %W %x %y"
	bind $widget <Any-ButtonRelease> "$this EndMotion %W %x %y"
	# SGS Don't allow rotation for 2D viewer
	#	bind $widget <B1-Motion> "$this Rotate %W %x %y"
	bind $widget <B2-Motion> "$this Pan %W %x %y"
	bind $widget <Shift-B1-Motion> "$this Pan %W %x %y"
	bind $widget <B3-Motion> "$this Zoom %W %x %y"
	
	bind $widget <KeyPress-c> {cvisTkCon}
	bind $widget <KeyPress-r> "$this Reset %W %x %y"
	
	bind $widget <KeyPress-w> "$this Wireframe %W"
	bind $widget <KeyPress-s> "$this Surface %W"
	
	bind $widget <KeyPress-p> "$this PickActor %W %x %y"
	
	set RenderWindow [$widget GetRenderWindow]

	$RenderWindow AddRenderer $Renderer

	set WindowNameProp [cvisProperty $name.WindowName]
	$name.WindowName AddCallback "$this UpdateWindowName"


	# SGS this was added to create a default light but 
	# there should be some better way of doing this?
	# If this is not here then we get no light and nothing
	# ever shows up
	$Renderer Render
	UpdateRenderer $RenderWidget 0 0 

	if { [file exists $::cvis::options::DEFAULT_CAMERA] } {
	    RestoreView	 $::cvis::options::DEFAULT_CAMERA
	}
    }

    destructor {

	$RenderWindow Delete
	$Renderer Delete
    }
    
    method GetPicker { } {
    }

    method GetRenderWindowInteractor { } {
    }

    method GetRenderer { } {
	return $Renderer
    }

    method GetWindowName {} {
	return [$WindowNameProp GetValue]
    }

    method GetWindowNameProp {} {
	return [$WindowNameProp GetReference]
    }

    method SetWindowName {windowname} {
	$WindowNameProp SetValue $windowname
    }
    
    method SetEquivWindowName { property } {
	$WindowNameProp SetEquiv $property
    }

    method UpdateWindowName { windowname } {
    }

    method GetRenderWindow { } {
	return $RenderWindow
    }
    
    method Render { } {
	if { [string length $CurrentLight] } {
	    # Make Light Follow Camera
	    eval $CurrentLight SetPosition [$CurrentCamera GetPosition]
	    eval $CurrentLight SetFocalPoint [$CurrentCamera GetFocalPoint]
	    $RenderWindow Render    
	}
    }

    method SaveImageAsPPM {} {
	# SGS
	# Warning don't check this in it will not work
	# VTK4 change needed
	# Use a vtkWindowToImageFilter instead of SaveImageAsPPM    
	#	$RenderWindow SaveImageAsPPM 
	# SGS
    }

    method SetFileName { filename } {
	$RenderWindow SetFileName $filename
    }

    method GetFileName { } {
	return [$RenderWindow GetFileName]
    }

    method SetSize { size_x size_y } {
	set Size(x) $size_x
	set Size(y) $size_y
	$RenderWindow SetSize $size_x $size_y
    }

    method SaveState { filename } {
	SaveView [$Renderer GetActiveCamera] $filename
    }

    method RestoreState { filename} {
	RestoreView [$Renderer GetActiveCamera] $filename
	Render
    }

    method AddActor { actor } {
	$Renderer AddActor $actor

    }

    method RemoveActor { actor } {
	$Renderer RemoveActor $actor

    }

    method AddVolume { volume } {
	$Renderer AddVolume $volume
    }

    method RemoveVolume { volume } {
	$Renderer RemoveVolume $volume
    }


    method SetUniformScaling {uniformscaling} {
	set UniformScaling $uniformscaling
	
	SetScaling
    }

    variable YLabel "Y"
    method SetYRange {min max delta label} {
	set YLabel $label
	
	$frame.yscroll configure \
		-from $min -to $max \
		-resolution $delta \
		-digits 9

	$frame.ylabel configure -text "$YLabel"
    }

    variable XLabel "X"
    method SetXRange {min max delta label} {
	set XLabel $label

	$frame.xscroll configure \
		-from $min -to $max \
		-resolution $delta \
		-digits 9

	$frame.xlabel configure -text "$XLabel"
    }

    variable ZLabel "Zoom"
    method SetZoomRange {min max delta base label} {
	set ZLabel $label

	set BaseZoom $base
	$frame.zoom configure \
		-from $min -to $max \
		-resolution $delta \
		-digits 9

	$frame.zlabel configure -text "$ZLabel"
    }

    method SetScaling { } {

	set collection [$Renderer GetActors]
	
	if {$UniformScaling } { 

	    # Reset scaling for bounds computation
	    # This does not seem to cause a speed problem but is
	    # there a better way to compute this
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		$item SetScale 1.0 1.0 1.0
		set item [$collection GetNextItem]
	    }
	    
	    # Scale display objects so they have the same size
	    # along each of their dimensions.  
	    set GlobalScaling $::cvis::math::FLT_MAX
	    foreach plane "x y z" { 
		set Delta($plane)  $::cvis::math::FLT_MAX
		set Bounds(min,$plane) $::cvis::math::FLT_MAX
		set Bounds(max,$plane) [expr -$::cvis::math::FLT_MAX]
	    }
	    
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {

		if { $GlobalScaling > [$item GetScaling] } {
		    set GlobalScaling [$item GetScaling]
		}
	    
		scan [$item GetBounds] "%f %f %f %f %f %f" \
			real(min,x) real(max,x) \
			real(min,y) real(max,y) \
			real(min,z) real(max,z)

		scan [$item GetSpacing] "%f %f %f" delta(x) delta(y) delta(z)
		
		foreach plane "x y z" {
		    if { $Delta($plane) > $delta($plane) } {
			set Delta($plane) $delta($plane)
		    }

		    if { $Bounds(min,$plane) > $real(min,$plane) } {
			set Bounds(min,$plane) $real(min,$plane)
		    }
		    
		    if { $Bounds(max,$plane) < $real(max,$plane) } {
			set Bounds(max,$plane) $real(max,$plane)
		    }
		}
		
		set item [$collection GetNextItem]
	    }

	    # SGS instead of scaling to 1000.0 what about scaling to
	    # Max extent size to avoid size changing as much

	    set Scaling(x) [expr 1000.0 / ( $Bounds(max,x) - $Bounds(min,x) )]
	    set Scaling(y) [expr 1000.0 / ( $Bounds(max,y) - $Bounds(min,y) )]
	    set Scaling(z) [expr 1000.0 / ( $Bounds(max,y) - $Bounds(min,y) )]
	} {
	    
	    # Do not do any scaling
	   
	    set Scaling(x) 1.0
	    set Scaling(y) 1.0
	    set Scaling(z) 1.0
	}

	# Set the scroll values

	# Convert back to original user coordinates
	foreach plane "x y z" {
	    set Delta($plane) [expr $Delta($plane) * $GlobalScaling]
	    set Bounds(min,$plane) [expr $Bounds(min,$plane) * $GlobalScaling]
	    set Bounds(max,$plane) [expr $Bounds(max,$plane) * $GlobalScaling]
	}



	# sgs fix this
	# What did I intend to go here?
	case $Plane { 
	    "z" {
		
	    }
	}
       
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
 	    $item SetScale $Scaling(x) $Scaling(y) $Scaling(z)
	    set item [$collection GetNextItem]
	}
    }

    variable AutoSave 0

    method SetAutoSave { autosave } {
	set AutoSave $autosave
    }

    method GetAutoSave { } {
	return $AutoSave
    }

    method Update {} {
	# SGS If autosave turned on then save the image to the current
	# file
	if {$AutoSave} { 
	    SaveImageAsPPM 
    }

    # Save the current state of the camera to a file.  
    # If not filename is passed in a file browser will be displayed
    method SaveView {{filename ""}} {

	set camera [[GetRenderer] GetActiveCamera]

	set types {
	    {{vtk Camera }  {.cam}        }
	    {{All Files }   *             }
	}
	
	if { $filename == "" } {
	    set filename [tk_getSaveFile -filetypes $types]
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
		puts $fp "Parallel Projection: \
			[expr int([$camera GetParallelProjection])]"
		puts $fp "Thickness: [$camera GetThickness]"
		puts $fp "View Angle: [$camera GetViewAngle]"
		puts $fp "View Plane Normal: [$camera GetViewPlaneNormal]"
		puts $fp "View Up: [$camera GetViewUp]"
		close $fp
	    }
	}
    }

    ## Recover view from saved file
    ## Provide camera name
    method RestoreView {{filename ""}} {

	set camera [[GetRenderer] GetActiveCamera]

	set types {
	    {{vtk Camera }  {.cam}        }
	    {{All Files }   *             }
	}
	
	if { $filename == "" } {
	    set filename [tk_getOpenFile -filetypes $types]
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
			# SGS Added this
			set YCenter [expr ($FPoint1 - $RPoint1)/2.0 + $FPoint1]
			set XCenter [expr ($FPoint0 - $RPoint0)/2.0 + $FPoint0]
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
		    
		    set loc [string first "View Plane Normal:" $line]
		    if { $loc != -1 } {
			set substr [string range $line [expr $loc + 18] end]
			scan $substr "%f %f %f" vpn1 vpn2 vpn3
			$camera SetViewPlaneNormal $vpn1 $vpn2 $vpn3
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

    method About { } {
	set tracervar [::itcl::local cvisTrace #auto]
	::cvis::gui::About
    }
}





