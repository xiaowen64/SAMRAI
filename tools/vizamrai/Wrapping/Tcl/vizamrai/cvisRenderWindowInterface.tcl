##
## File:        cvisRenderWindowInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Interface for Render Window for 3D data 
##

class cvisRenderWindowInterface {
    inherit cvisInterface

    variable FileName ""
    variable FileNameProp

    variable CameraFileName ""
    variable CameraFileTypes {{Camera {.cam}}}

    # SGS if the scaling is true then loop over the actors and 
    # scale by the correct amount?
    variable UniformScaling 0
    variable ParallelProjection 0
    variable Size
    variable OffScreenRendering 0

    variable AutoSave 0

    variable Window ""

    variable PickPointIsVisible 1

    variable Lights 

    constructor {name} {
	cvisInterface::constructor $name
    } {
	set FileNameProp [cvisProperty $name.FileName]
	$name.FileName AddCallback "$this UpdateFileName"

	set Size(x) 640
	set Size(y) 480

	foreach light {1 2 3 4} {
	    set Lights($light,intensity) 0.0
	    set Lights($light,x) 0.0
	    set Lights($light,y) 0.0
	    set Lights($light,z) 0.0
	    set Lights($light,vtklight) 0
	}
	
	set light camera
	set Lights($light,intensity) 1.0
	set Lights($light,x) 0.0
	set Lights($light,y) 0.0
	set Lights($light,z) 0.0
    }

    method Save {file} {
	if {[string length $FileName]} {
	    puts $file "$this SetFileName \"$FileName\""
	}
	puts $file "$this SetUniformScaling $UniformScaling"
	puts $file "$this SetParallelProjection $ParallelProjection"
	puts $file "$this SetAutoResize $AutoResize"
	puts $file "$this SetSize $Size(x) $Size(y)"
	puts $file "$this SetOffScreenRendering $OffScreenRendering"
	puts $file "$this SetAutoSave $AutoSave"
	puts $file "$this SetBackgroundColor $BackgroundColor"
	if {[string length $CameraFileName]} {
	    puts $file "$this SetCameraFileName \"$CameraFileName\""
	    $Window SaveView $CameraFileName
	}
	puts $file "$this SetPickPointVisibility $PickPointIsVisible"
	puts $file "$this SetPickPointColor $PickPointColor"
    }

    method ExecuteSetLight {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {[string length $Window]} {
	    
	    # First do special stuff with the camera light
	    set camera_light [$Window GetLight camera]
	    $camera_light SetIntensity $Lights(camera,intensity)

	    # Process the other lights in the scene
	    foreach light {1 2 3 4} {
		set light_obj [$Window GetLight $light]
		if { $Lights($light,intensity) > 0.0 } {
		    $light_obj SetSwitch 1
		    $light_obj SetIntensity $Lights($light,intensity)
		    $light_obj SetPosition \
			    $Lights($light,x) \
			    $Lights($light,y) \
			    $Lights($light,z)
		} {
		    $light_obj SetSwitch 0
		}
	    }

	    $Window Render
	}
    }
    
    method SetRenderWindow {window} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Window $window
	$window SetInterface $this
    }

    method GetRenderWindow {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Window
    }

    method GetFileName { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$FileNameProp GetValue]
    }
    
    method GetFileNameProp {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$FileNameProp GetReference]
    }

    variable CameraFileBrowser ".null"

    method SetCameraFileName {camerafilename} {
	set CameraFileName $camerafilename
	$Window RestoreView $CameraFileName
    }

    method SaveCamera {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $CameraFileBrowser]]} {
	    $CameraFileBrowser SetSavePrompt
	} {
	    set CameraFileBrowser [cvisFileBrowser #auto]
	    $CameraFileBrowser SetTitle "Select Camera File"
	    $CameraFileBrowser SetFileTypes $CameraFileTypes
	    $CameraFileBrowser SetSavePrompt
	}

	set CameraFileName [$CameraFileBrowser Prompt]

	$Window SaveView $CameraFileName
    }

    method RestoreCamera {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $CameraFileBrowser]]} {
	    $CameraFileBrowser SetOpenPrompt
	} {
	    set CameraFileBrowser [cvisFileBrowser #auto]
	    $CameraFileBrowser SetTitle "Select Camera File"
	    $CameraFileBrowser SetFileTypes $CameraFileTypes
	    $CameraFileBrowser SetOpenPrompt
	}

	set CameraFileName [$CameraFileBrowser Prompt]

	$Window RestoreView $CameraFileName
    }

    method SetFileName { filename } {
	set tracervar [::itcl::local cvisTrace #auto]
	set FileName $filename
	$FileNameProp SetValue $filename
    }

    method SetEquivWindowName { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$FileNameProp SetEquiv $property
    }

    method UpdateFileName {filename} {
	set tracervar [::itcl::local cvisTrace #auto]
	# SGS why did we do this?
	set temp "[file rootname $filename].[[$Window GetImageWriter] GetType]"
	set FileName $temp
	$Window SetFileName $FileName
    }

    variable ImageFileBrowser ".null"

    method FileNameBrowser { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length [find objects $ImageFileBrowser]]} {
	    $ImageFileBrowser SetSavePrompt
	} {
	    set ImageFileBrowser [cvisFileBrowser #auto]
	    $ImageFileBrowser SetTitle "Select Image File"
	    $ImageFileBrowser SetFileTypes [[$Window GetImageWriter] GetFileTypes]
	    $ImageFileBrowser SetSavePrompt
	}

	$ImageFileBrowser SetDefaultExtension [[$Window GetImageWriter] GetType]

	set FileName [$ImageFileBrowser Prompt]

	if {[string length $FileName]} {
	    # Set the file type based on the filename the user selected
	    [$Window GetImageWriter] SetType [string trimleft [file extension $FileName] "."]
	    SetFileName $FileName
	}
    }

    method SetUniformScaling {uniformscaling } {
	if {$UniformScaling != $uniformscaling } {
	    set UniformScaling $uniformscaling
	    ExecuteSetUniformScaling
	}
    }

    method GetUniformScaling {} {
	return $UniformScaling
    }

    method ExecuteSetUniformScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Window] } { 
	    $Window SetUniformScaling $UniformScaling
	    $Window Render
	}
    }

    method SetParallelProjection {parallelprojection } {
	if {$ParallelProjection != $parallelprojection } {
	    set ParallelProjection $parallelprojection
	    ExecuteSetParallelProjection
	}
    }

    method GetParallelProjection {} {
	return $ParallelProjection
    }

    method ExecuteSetParallelProjection {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Window] } { 
	    $Window SetParallelProjection $ParallelProjection
	    $Window Render
	}
    }

    variable AutoResize 1
    method SetAutoResize {resize} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$AutoResize != $resize } {
	    set AutoResize $resize
	    ExecuteSetAutoResize
	}
    }

    method GetAutoResize {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $AutoResize
    }

    method ExecuteSetAutoResize {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Window] } { 
	    $Window SetAutoResize $AutoResize
	}
    }
    

    method SetSize {x y} {
	set Size(x) $x
	set Size(y) $y
	ExecuteSetSize
    }

    method GetSize {} {
	return $Size(x) $Size(y)
    }

    method ExecuteSetSize {} {
	if {[string length $Window]} {
	    $Window SetSize $Size(x) $Size(y)
	}
    }

    method SetOffScreenRendering {offscreenrendering} {
	if {$OffScreenRendering != $offscreenrendering} {
	    set OffScreenRendering $offscreenrendering
	    ExecuteSetOffScreenRendering
	}
    }

    method ExecuteSetOffScreenRendering {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Window] } { 
	    $Window SetOffScreenRendering $OffScreenRendering
	}
    }


    method SetAutoSave {autosave} {
	if {$AutoSave != $autosave} {
	    set AutoSave $autosave
	    ExecuteSetAutoSave
	}
    }

    method ExecuteSetAutoSave {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Window] } { 
	    $Window SetAutoSave $AutoSave
	}
    }

    method Write { } {
	set tracervar [::itcl::local cvisTrace #auto]
	$Window Write
    }

    variable BackgroundColor "0.0 0.0 0.0"

    method SetBackgroundColor {r g b} {
	set BackgroundColor "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	
	set frame $NotebookPage.optionsframe
	set cs [$frame childsite]
	set frame $cs.backgroundColor
	$frame.color configure -bg $color
	
	$Window SetBackground $r $g $b
    }

    method SelectBackgroundColor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $BackgroundColor
	
	set value [color Prompt]
	
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetBackgroundColor $r $g $b
	}	

	delete object color 
    }


    variable PickPointColor "1.0 1.0 1.0"

    method SetPickPointColor {r g b} {
	set PickPointColor "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	
	set frame $NotebookPage.pickpoint
	set cs [$frame childsite]
	$cs.color.color configure -bg $color

	if {[string length $Window]} {
	    $Window SetPickPointColor $r $g $b
	    $Window Render
	}
    }

    method SelectPickPointColor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $PickPointColor
	
	set value [color Prompt]
	
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetPickPointColor $r $g $b
	}	

	delete object color 
    }

    method ExecuteSetPickPointVisibility {} {
	if {[string length $Window]} {
	    $Window SetPickPointVisibility $PickPointIsVisible
	    $Window Render
	}
    }

    method SetPickPointVisibility {visible} {
	if {$visible != $PickPointIsVisible} {
	    set PickPointIsVisible $visible
	}
    }
    
    method GetPickPointVisibility {} {
	return $PickPointIsVisible
    }

    method SetNotebookPage {notebook} {
	set tracervar [::itcl::local cvisTrace #auto]
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	# Construct controls for image output
	#-----------------------------------------------------------------
	set frame $NotebookPage.filenameframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Image Output"
	set cs [$frame childsite]

	button $cs.save \
		-text "Save Current Image" \
		-command "$this Write"
	
	button $cs.selectfile \
		-text "FileName" \
		-command  "$this FileNameBrowser" 

	entry $cs.filename -textvariable [scope FileName]

	checkbutton $cs.autosave \
		-text "Automatically Save on New Load" -variable [scope AutoSave] \
		-command "$this ExecuteSetAutoSave"

	checkbutton $cs.offscreenrendering \
		-text "Off Screen Rendering" \
		-variable [scope OffScreenRendering] \
		-command "$this ExecuteSetOffScreenRendering"
	
	#-----------------------------------------------------------------
	# Construct controls for render options
	#-----------------------------------------------------------------
	set frame $NotebookPage.optionsframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Render Options"
	set cs [$frame childsite]

	checkbutton $cs.uniformscale \
		-text "Uniformly Scale Axis" -variable [scope UniformScaling] \
		-command "$this ExecuteSetUniformScaling"

	checkbutton $cs.parallelprojection \
		-text "Parallel Projection" -variable [scope ParallelProjection] \
		-command "$this ExecuteSetParallelProjection"

	checkbutton $cs.autoresize \
		-text "Grow render area with window" -variable [scope AutoResize] \
		-command "$this ExecuteSetAutoResize"

	set frame $cs.size
	frame $frame

	label $frame.widthtitle -text "Width"
	entry $frame.width \
		-textvariable [scope Size(x)] 
	bind $frame.width <Return> "$this ExecuteSetSize"

	label $frame.heighttitle -text "Height"
	entry $frame.height \
		-textvariable [scope Size(y)] 
	bind $frame.height <Return> "$this ExecuteSetSize"

	set frame $cs.backgroundColor
	frame $frame
	label $frame.label -text "Background Color"
	
	scan $BackgroundColor "%f %f %f" r g b
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	button $frame.color -bg $color \
		-command "$this SelectBackgroundColor"

	#-----------------------------------------------------------------
	# Save/Restore camera view
	#-----------------------------------------------------------------
	set frame $NotebookPage.saverestore
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Save and Restore Camera View"
	set cs [$frame childsite]

	button $cs.savecamera -text "Save" \
		-command  "$this SaveCamera" 

	button $cs.restorecamera -text "Restore" \
		-command  "$this RestoreCamera" 

	entry $cs.filename \
		-textvariable [scope CameraFileName]

	#-----------------------------------------------------------------
	#
	#-----------------------------------------------------------------
	set frame $NotebookPage.pickpoint
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Pick Point"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "PickPoint Visible" -variable [scope PickPointIsVisible] \
		-command "$this ExecuteSetPickPointVisibility"

	frame $cs.color

	label $cs.color.label -text PickPointColor
	
	scan $PickPointColor "%f %f %f" r g b
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	button $cs.color.color -bg $color \
		-command "$this SelectPickPointColor"

	#-----------------------------------------------------------------
	# Lights
	#-----------------------------------------------------------------
	set frame $NotebookPage.lights
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Lighting"
	set cs [$frame childsite]

	set lframe $cs.lights
	frame $lframe
	set light "title"
	label $lframe.name_$light -text "Name"
	label $lframe.intensity_$light -text "Intensity"
	label $lframe.x_$light -text "X"
	label $lframe.y_$light -text "Y"
	label $lframe.z_$light -text "Z"
	foreach light {camera 1 2 3 4} {
	    label $lframe.name_$light -text "$light"
	    scale $lframe.intensity_$light \
		    -from 0 \
		    -to 1 \
		    -resolution 0.05 \
		    -variable [scope Lights($light,intensity)] \
		    -orient horizontal \
		    -highlightthickness 0
	    bind  $lframe.intensity_$light <ButtonRelease> "$this ExecuteSetLight"

	    if {[string compare $light "camera"]} {
		scale $lframe.x_$light \
			-from -1 \
			-to 1 \
			-resolution 0.1 \
			-variable [scope Lights($light,x)] \
			-orient horizontal \
			-highlightthickness 0
		bind  $lframe.x_$light <ButtonRelease> "$this ExecuteSetLight"
		
		scale $lframe.y_$light \
			-from -1 \
			-to 1 \
			-resolution 0.1 \
			-variable [scope Lights($light,y)] \
			-orient horizontal \
			-highlightthickness 0
		bind  $lframe.y_$light <ButtonRelease> "$this ExecuteSetLight"
		
		scale $lframe.z_$light \
			-from -1 \
			-to 1 \
			-resolution 0.1 \
			-variable [scope Lights($light,z)] \
			-orient horizontal \
		    -highlightthickness 0
		bind  $lframe.z_$light <ButtonRelease> "$this ExecuteSetLight"
	    }
	}

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.filenameframe
	set cs [$frame childsite]
	pack $cs.save -side top -anchor w \
		-padx 4 -pady 8
	pack $cs.autosave  -anchor w -pady 4
	pack $cs.offscreenrendering  -anchor w -pady 4
	pack $cs.selectfile -side left \
		-padx 4 -pady 4
	pack $cs.filename -side left \
		-fill x -expand yes -pady 4
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.optionsframe
	pack $frame -fill x
	set cs [$frame childsite]
	pack $cs.uniformscale -anchor w -pady 4
	pack $cs.parallelprojection -anchor w
	pack $cs.autoresize -anchor w

	set frame $cs.size
	pack $frame.widthtitle -side left 
	pack $frame.width -side left -padx 4 -pady 4
	pack $frame.heighttitle -side left
	pack $frame.height -side left -padx 4 -pady 4
	pack $frame -fill x

	set frame $cs.backgroundColor
	pack $frame.label -side left
	pack $frame.color -side right
	pack $frame -anchor w

	pack $cs -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.saverestore
	set cs [$frame childsite]
	pack $cs.savecamera -side left  -anchor w -padx 4 -pady 8
	pack $cs.restorecamera -side left -anchor w -padx 4 -pady 4
	pack $cs.filename -side left \
		-fill x -expand yes -pady 4
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.pickpoint
	set cs [$frame childsite]
	pack $cs.visible -anchor w -pady 4
	pack $cs.color.label -side left
	pack $cs.color.color -side right
	pack $cs.color -anchor w
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.lights
	set cs [$frame childsite]

	set lframe $cs.lights
	foreach light {title camera 1 2 3 4} row {0 1 2 3 4 5} {
	    grid $lframe.name_$light -row $row -column 0 -pady 4
	    grid $lframe.intensity_$light -row $row -column 1
	    if {[string compare $light "camera"]} {
		grid $lframe.x_$light -row $row -column 2
		grid $lframe.y_$light -row $row -column 3
		grid $lframe.z_$light -row $row -column 4
	    }
	}
	pack $lframe -fill x

	pack $cs -fill x
	pack  $frame -fill x 
    }
}







