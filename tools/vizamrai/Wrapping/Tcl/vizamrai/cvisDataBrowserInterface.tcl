##
## File:        cvisDataBrowserInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: "Real" Data Browser Interface
##


class cvisDataBrowserInterface {
    inherit cvisInterface

    variable RenderWindowIsVisible 0

    variable DataIsVisible 1
    variable LabelIsVisible 1

    variable CellBoundaryWidth 3
    variable CellBoundaryIsVisible 1

    variable Input ""

    variable DataLabelFormat "%1.2e"
    variable AxisLabelFormat "%1.2e"

    variable Slice 0
    variable SlicePlane "z"
    variable PickPoint
    variable PickPointReal
    variable WindowSize

    variable RealBounds 
    variable IndexBounds
    variable IndexScaling
    variable Scaling

    constructor {name} {
	cvisInterface::constructor $name
    } {
	foreach plane {x y z} {
	    set RealBounds(min,$plane) 0.0
	    set RealBounds(max,$plane) 1.0

	    set IndexBounds(min,$plane) 0
	    set IndexBounds(max,$plane) 1

	    set IndexScaling($plane) 1.0
	    set PickPoint($plane) 1
	    set PickPointReal($plane) 0.5
	    set WindowSize($plane) 4
	}

    }

    method Save {file} {
	set tracervar [::itcl::local cvisTrace #auto]
	puts $file "$this SetDataFormat \"[GetDataFormat]\""
	puts $file "$this SetAxisFormat \"[GetAxisFormat]\""
	puts $file "$this SetRenderWindowVisibility  [GetRenderWindowVisibility]"
	puts $file "$this SetWindowSize [GetWindowSize]"
	puts $file "$this SetPickPoint  [GetPickPoint]"
	puts $file "$this SetSlicePlane [GetSlicePlane]"
    }


    method ExecuteSetDataFormat {} {
	foreach plane $Input {
	    $plane SetDataLabelFormat $DataLabelFormat
	}
    }
    
    method SetDataFormat {format} {
	set DataLabelFormat $format
	ExecuteSetDataFormat
    }

    method GetDataFormat {} {
	return $DataLabelFormat
    }

    method GetPickPoint {} {
	return "$PickPointReal(x) $PickPointReal(y) $PickPointReal(z)"
    }

    # input is in screen coordinates
    method SetPickPoint {x y z} {

	set PickPointReal(x) [expr $x / $Scaling]
	set PickPointReal(y) [expr $y / $Scaling]
	set PickPointReal(z) [expr $z / $Scaling]

	set frame $NotebookPage.sliceframe
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set PickPoint($plane) \
		    [expr int( ($PickPointReal($plane) - $RealBounds(min,$plane)) / $IndexScaling($plane))]
	    set PickPointReal($plane) [expr ($PickPoint($plane)+0.5) * \
		    $IndexScaling($plane) + $RealBounds(min,$plane)]
	    set screen($plane) [expr $PickPointReal($plane) * $Scaling]
	    $cs.$plane.realvalue configure \
		    -text [format "%1.3e" $PickPointReal($plane)]
	}

	foreach item $Input {
	    $item SetPickPoint $screen(x) $screen(y) $screen(z)
	}
    }
    
    method ExecuteSetWindowSize {} {
	set tracervar [::itcl::local cvisTrace #auto]	
	foreach item $Input {
	    $item SetWindowSize $WindowSize(x) $WindowSize(y) $WindowSize(z)
	}

    }

    method SetWindowSize {x y z} {
	set WindowSize(x) $x
	set WindowSize(y) $y
	set WindowSize(z) $z
    }

    method GetWindowSize {} {
	return "$WindowSize(x) $WindowSize(y) $WindowSize(z)"
    }

    method PickPointChange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set frame $NotebookPage.sliceframe
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set real($plane) [expr ($PickPoint($plane)+0.5) * \
		    $IndexScaling($plane) + $RealBounds(min,$plane)]
	    set screen($plane) [expr $real($plane) * $Scaling]
	    $cs.$plane.realvalue configure \
		    -text [format "%1.3e" $real($plane)]
	}

	foreach item $Input {
	    $item SetPickPoint $screen(x) $screen(y) $screen(z)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow SetPickPoint $screen(x) $screen(y) $screen(z)
	}
    }

    method ExecuteSetSlicePlane {} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach item $Input {
	    $item SetSlicePlane $SlicePlane
	}
    }

    method SetSlicePlane {plane} {
	if {$plane != $SlicePlane} {
	    set SlicePlane $plane
	    ExecuteSetSlicePlane 
	}
    }

    method GetSlicePlane {} {
	return $SlicePlane
    }

    method ExecuteSetAxisFormat {} {
	foreach plane $Input {
	    $plane SetAxisLabelFormat $AxisLabelFormat
	}
    }

    method SetAxisFormat {format} {
	set AxisLabelFormat $format
	ExecuteSetAxisFormat
    }

    method GetAxisFormat {} {
	return $AxisLabelFormat
    }

    method AddInput {browser} {
	lappend Input $browser

	$browser SetInterface $this
	
	Update
    }

    method ExecuteSetRenderWindowVisibility {} {
	foreach plane $Input {
	    $plane SetRenderWindowVisibility $RenderWindowIsVisible
	}
    }

    method SetRenderWindowVisibility {visible} {
	if {$visible != $RenderWindowIsVisible} {
	    set RenderWindowIsVisible $visible
	    ExecuteSetRenderWindowVisibility
	}
    }

    method GetRenderWindowVisibility {} {
	return $RenderWindowIsVisible
    }

    method Update {} {

	# SGS what do we do if multiple inputs?
	set item [lindex $Input 0]

	set Scaling [$item GetScaling]

	scan [$item GetIndexBounds] "%f %f %f %f %f %f" \
		IndexBounds(min,x) IndexBounds(max,x) \
		IndexBounds(min,y) IndexBounds(max,y) \
		IndexBounds(min,z) IndexBounds(max,z)

	scan [$item GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z) 

	scan [$item GetBounds] "%f %f %f %f %f %f" \
		RealBounds(min,x) RealBounds(max,x) \
		RealBounds(min,y) RealBounds(max,y) \
		RealBounds(min,z) RealBounds(max,z)

	# convert from screen to real coords
	foreach plane {x y z} {
	    set RealBounds(min,$plane) [expr $RealBounds(min,$plane) / $Scaling]
	    set RealBounds(max,$plane) [expr $RealBounds(max,$plane) / $Scaling]
	}

	set frame $NotebookPage.sliceframe
	set cs [$frame childsite]

	foreach plane {x y z} {
	    $cs.$plane.scale configure \
		    -from $IndexBounds(min,$plane) -to $IndexBounds(max,$plane)
	}
    }

    method SetNotebookPage {notebook} {
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	#
	#-----------------------------------------------------------------
	set frame $NotebookPage.window
	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Data Browser Window"
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken

	#-----------------------------------------------------------------
	# Visibility
	#-----------------------------------------------------------------
	set frame $NotebookPage.visible

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Window Visibility"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "Window Visible" -variable [scope RenderWindowIsVisible] \
		-command "$this ExecuteSetRenderWindowVisibility"

	#-----------------------------------------------------------------
	# X/Y/Z Plane Selector
	#-----------------------------------------------------------------
	set frame $NotebookPage.planeframe

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Slice Plane"
	set cs [$frame childsite]

	frame $cs.plane

	radiobutton $cs.plane.x -text "X" \
		-variable [scope SlicePlane] \
		-value "x" \
		-command "$this ExecuteSetSlicePlane"

	radiobutton $cs.plane.y -text "Y" \
		-variable [scope SlicePlane] \
		-value "y" \
		-command "$this ExecuteSetSlicePlane"

	radiobutton $cs.plane.z -text "Z" \
		-variable [scope SlicePlane] \
		-value "z" \
		-command "$this ExecuteSetSlicePlane"

	#-----------------------------------------------------------------
	# PointSelector
	#-----------------------------------------------------------------
	set frame $NotebookPage.sliceframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext  "Center Point"
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set frame $cs.$plane
	    frame $frame
	    label $frame.slicename \
		    -text "$plane"

	    set real($plane) [expr ($PickPoint($plane)+0.5) * $IndexScaling($plane) + $RealBounds(min,$plane)]
	    label $frame.realvalue \
		    -text [format "%1.3e" $real($plane)]

	    scale $frame.scale \
		    -length 100 \
		    -from $IndexBounds(min,$plane) \
		    -to $IndexBounds(max,$plane) \
		    -variable [scope PickPoint($plane)] \
		    -orient horizontal \
		    -highlightthickness 0
	    
	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $frame.scale <ButtonRelease> \
		    "$this PickPointChange"
	}

	#-----------------------------------------------------------------
	# Window Size
	#-----------------------------------------------------------------
	set frame $NotebookPage.sizeframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Window Size"
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set frame $cs.$plane
	    frame $frame
	    label $frame.slicename \
		    -text "$plane"

	    scale $frame.scale \
		    -from 1 \
		    -to 20 \
		    -variable [scope WindowSize($plane)] \
		    -orient horizontal \
		    -highlightthickness 0
	    
	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $frame.scale <ButtonRelease> \
		    "$this ExecuteSetWindowSize"
	}

	#-----------------------------------------------------------------
	# Text Label
	#-----------------------------------------------------------------
	set frame $NotebookPage.label
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Labels"
	set cs [$frame childsite]


	set frame $cs.data
	frame $frame 

	label $frame.labelformattext \
		-text "Data Label Format"
	
	entry $frame.labelformat \
		-textvariable [scope DataLabelFormat]

	bind $frame.labelformat <Return> "$this ExecuteSetDataFormat"

	set frame $cs.axis
	frame $frame

	label $frame.labelformattext \
		-text "Axis Label Format"
	
	entry $frame.labelformat \
		-textvariable [scope AxisLabelFormat]

	bind $frame.labelformat <Return> "$this ExecuteSetAxisFormat"

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.window
	pack $frame.title
	pack $frame.titleseperator -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.visible
	set cs [$frame childsite]

	pack $cs.visible -anchor w -pady 4
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.planeframe
	set cs [$frame childsite]

	pack $cs.plane -expand yes -fill x 
	pack $cs.plane.x -side left -pady 4
	pack $cs.plane.y -side left
	pack $cs.plane.z -side left

  	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.sliceframe
	pack $frame -fill x
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set frame $cs.$plane
	    pack $frame.slicename -side left
	    pack $frame.realvalue -side right 
	    pack $frame.scale -expand yes -fill x -pady 2
	    
	    pack $frame -fill x
	}

	pack $cs -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.sizeframe
	pack $frame -fill x
	set cs [$frame childsite]

	foreach plane {x y z} {
	    set frame $cs.$plane
	    pack $frame.slicename -side left
	    pack $frame.scale -expand yes -fill x 
	    
	    pack $frame -fill x
	}

	#-----------------------------------------------------------------
	set frame $NotebookPage.label
	set cs [$frame childsite]

	pack $cs.data.labelformattext -anchor w -side left
	pack $cs.data.labelformat -anchor w 
	pack $cs.data -fill x
	
	pack $cs.axis.labelformattext -anchor w  -side left
	pack $cs.axis.labelformat -anchor w 
	pack $cs.axis -fill x

	pack $frame -fill x
    }   
}

