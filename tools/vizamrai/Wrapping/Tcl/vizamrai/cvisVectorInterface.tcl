##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisVectorInterface.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Slice plane interface
##

class cvisVectorInterface {
    inherit cvisInterface

    common Number 0
    
    variable VOI "0 1 0 1 0 1"
    variable Corner
    variable Depth
    
    variable Crop 0
    variable CroppingRange 

    variable Input ""
    variable PlaneInput ""

    variable Scaling 1.0

    variable IndexScaling
    variable MinIndexScaling
    variable IndexBounds

    variable RealBounds
    
    variable TableRange

    variable DataIsVisible 1

    variable ColorMap 1
    variable VectorColor "1 1 1"

    variable VectorScaleFactor

    method Save {file} {
	puts $file "$this SetColorMap [GetColorMap]"
	puts $file "$this SetVectorColor [GetVectorColor]"
	puts $file "$this SetVOI [GetVOI]"
	puts $file "$this SetCrop [GetCrop]"
	puts $file "$this SetCroppingRange [GetCroppingRange]"
	puts $file "$this SetDataVisible $DataIsVisible"
	puts $file "$this SetVectorScaleFactor [GetVectorScaleFactor]"

    }

    constructor {name} {
	cvisInterface::constructor $name
    } {

	foreach plane {x y z} {
	    set IndexBounds(min,$plane) 0
	    set IndexBounds(max,$plane) 0
	    set IndexScaling($plane) 1.0

	    set RealBounds(min,$plane) 0
	    set RealBounds(max,$plane) 1
	}

	set Corner(x) 0
	set Corner(y) 0
	set Corner(z) 0

	set Depth(x) 1
	set Depth(y) 1
	set Depth(z) 1

	set CroppingRange(min) 0.0
	set CroppingRange(max) 1.0

	set TableRange(min) 0.0
	set TableRange(max) 1.1

	set VectorScaleFactor(min) 0
	set VectorScaleFactor(max) 1
    }

    method ColorMapChange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane $Input {
	    $plane SetScalarVisibility $ColorMap
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }


    method GetColorMap {} {
	return $ColorMap
    }

    method SetColorMap {colormap} {
	set tracervar [::itcl::local cvisTrace #auto]
	set ColorMap $colormap
	ColorMapChange
    }

    method GetVectorColor {} {
	return $VectorColor
    }

    method SetVectorColor {r g b} {
	set VectorColor "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]

	set frame $NotebookPage.colorframe
	set cs [$frame childsite]
	$cs.vectorcolor configure -bg $color
	
	foreach plane $Input {
	    eval $plane SetVectorColor $r $g $b
	}
	
	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method ExecuteSetVectorColor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $VectorColor
	
	set value [color Prompt]
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetVectorColor $r $g $b
	}
	
	delete object color 
    }


    method AddPlaneInput {plane} {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend PlaneInput $plane
    }

    method AddInput {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]

	# SGS what to set VOI to?
	# set Slice 0

	lappend Input $sliceplane

	Update

	$sliceplane SetInterface $this

	VOIChange
    }

    method SetDataVisible {visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	set DataIsVisible $visible
	DataVisibleChange
    }

    method DataVisibleChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetDataVisibility $DataIsVisible
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    #-----------------------------------------------------------------------------
    method GetVOI {} {
	return $VOI
    }

    method SetVOI {min_x max_x min_y max_y min_z maz_z} {
	set tracervar [::itcl::local cvisTrace #auto]
	set VOI "$min_x $max_x $min_y $max_y $min_z $maz_z"

	VOIChange
    }

    method VOIChange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane {x y z} {
	    set real_corner($plane) \
		    [expr (($Corner($plane)+0.5) * $IndexScaling($plane)) \
		    + $RealBounds(min,$plane) ]
	    set screen_corner($plane) [expr $real_corner($plane) * $Scaling]

	    set real_depth($plane) \
		    [expr (($Depth($plane)-1) * $IndexScaling($plane))]

	    set screen_depth($plane) \
		    [expr $real_depth($plane) * $Scaling]

	    set frame $NotebookPage.corner
	    set cs [$frame childsite]
	    $cs.$plane.realvalue configure -text $real_corner($plane)

	    set frame $NotebookPage.depth
	    set cs [$frame childsite]
	    $cs.$plane.realvalue configure -text $real_depth($plane)
	}

	foreach plane $Input {
	    $plane SetVOI \
		    $screen_corner(x) [expr $screen_corner(x) + $screen_depth(x)] \
		    $screen_corner(y) [expr $screen_corner(y) + $screen_depth(y)] \
		    $screen_corner(z) [expr $screen_corner(z) + $screen_depth(z)]
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    #-----------------------------------------------------------------------------
    method GetCrop {} {
	return $Crop
    }

    method SetCrop {flag} {
	set Crop $flag
	CropChange
    }
    
    method CropChange {} {
	set item [lindex $Input 0]
	$item SetCrop $Crop

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    #-----------------------------------------------------------------------------
    method GetCroppingRange {} {
	return "$CroppingRange(min) $CroppingRange(max)"
    }

    method SetCroppingRange {min max} {
	set CroppingRange(min) $min
	set CroppingRange(max) $max
	CroppingRangeChange
    }

    method CroppingRangeChange {} {
	set item [lindex $Input 0]

	if {$CroppingRange(min) > $CroppingRange(max)} {
	    set CroppingRange(min) $CroppingRange(max)
	}

	$item SetCroppingRange $CroppingRange(min) $CroppingRange(max)

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    #-----------------------------------------------------------------------------
    method GetVectorScaleFactor {} {
	return "$VectorScaleFactor(min) $VectorScaleFactor(max)"
    }

    method SetVectorScaleFactor {min max} {
	set VectorScaleFactor(min) $min
	set VectorScaleFactor(max) $max
	VectorScaleFactorChange
    }

    method VectorScaleFactorChange {} {
	set screen_min [expr $VectorScaleFactor(min) * $MinIndexScaling * $Scaling]
	set screen_max [expr $VectorScaleFactor(max) * $MinIndexScaling *$Scaling]

	set item [lindex $Input 0]	
	$item SetVectorScaleFactor $screen_min $screen_max

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    #-----------------------------------------------------------------------------
    method GetSliceInScreenCoord {} {
	set real [expr (($Slice) * $IndexScaling($SlicePlane)) \
		+ $RealBounds(min,$SlicePlane) ]
	set screen [expr $real * $Scaling]
	return $screen
    }


    # SGS should update interface only, not planes pointed to
    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set item [lindex $Input 0]

	set Scaling [$item GetScaling]

	array set OldIndexBounds [array get IndexBounds]

	scan [$item GetIndexBounds] "%f %f %f %f %f %f" \
		IndexBounds(min,x) IndexBounds(max,x) \
		IndexBounds(min,y) IndexBounds(max,y) \
		IndexBounds(min,z) IndexBounds(max,z)

	foreach plane {x y z} {
	}

	scan [$item GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z)

	set MinIndexScaling $::cvis::math::FLT_MAX
	foreach plane {x y z} {
	    if {$MinIndexScaling > $IndexScaling($plane) } {
		set MinIndexScaling $IndexScaling($plane)
	    }
	} 

	scan [$item GetBounds] "%f %f %f %f %f %f" \
		RealBounds(min,x) RealBounds(max,x) \
		RealBounds(min,y) RealBounds(max,y) \
		RealBounds(min,z) RealBounds(max,z)

	scan [$item GetTableRange] "%f %f" \
		TableRange(min) TableRange(max)

	foreach plane {x y z} {
	    # convert from screen to real coords
	    set RealBounds(min,$plane) [expr $RealBounds(min,$plane) / $Scaling]
	    set RealBounds(max,$plane) [expr $RealBounds(max,$plane) / $Scaling]
	    
	    # Selection Box widgets
	    set frame $NotebookPage.corner
	    set cs [$frame childsite]
	    $cs.$plane.scale configure \
		    -from $IndexBounds(min,$plane) \
		    -to   $IndexBounds(max,$plane)

	    if {$OldIndexBounds(min,$plane) == $OldIndexBounds(max,$plane)} {
		set Corner($plane) $IndexBounds(min,$plane)
	    } {
		set Corner($plane) [expr int($Corner($plane) * \
			(($IndexBounds(max,$plane) - $IndexBounds(min,$plane)) / \
			($OldIndexBounds(max,$plane) - $OldIndexBounds(min,$plane)))) + \
			$OldIndexBounds(min,$plane)]
	    }

	    set old_depth $Depth($plane)

	    set frame $NotebookPage.depth
	    set cs [$frame childsite]
	    $cs.$plane.scale configure \
		    -from 1  \
		    -to [expr $IndexBounds(max,$plane) - $IndexBounds(min,$plane) + 1]

	    if {$OldIndexBounds(min,$plane) == $OldIndexBounds(max,$plane)} {
		set Depth($plane) [expr $IndexBounds(max,$plane) - $IndexBounds(min,$plane) + 1]
	    } {
		set Depth($plane) [expr int($old_depth * \
			(($IndexBounds(max,$plane) - $IndexBounds(min,$plane)) / \
			($OldIndexBounds(max,$plane) - $OldIndexBounds(min,$plane))))]
	    }
	}

	# Cropping widgets
	set frame $NotebookPage.cropping
	set cs [$frame childsite]
	foreach plane {min max} {
	    $cs.$plane.scale configure \
		    -from $TableRange(min) \
		    -to $TableRange(max) \
		    -resolution [expr ($TableRange(max) - $TableRange(min)) / 100]
	}

	# Vector Scale Factor widgets
	VectorScaleFactorChange
	VOIChange
    }
    

    # SGS need to have methods to set the scaler precision etc

    # SGS method to set the min and max scalar or always do the
    # add input method?

    method SetNotebookPage {notebook} {
	set tracervar [::itcl::local cvisTrace #auto]
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	# Slice Plane controls
	#-----------------------------------------------------------------

	set frame $NotebookPage.slice
	frame $frame -relief sunken -borderwidth 2

	label $frame.title -text "Slice Plane"
	frame $frame.titleseperator \
		-height 2 -borderwidth 1 -relief sunken

	#-----------------------------------------------------------------
	# Visibility
	#-----------------------------------------------------------------

	set frame $NotebookPage.visibility
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Visibility"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "Visible" -variable [scope DataIsVisible] \
		-command "$this DataVisibleChange"

	#-----------------------------------------------------------------
	# Corner
	#-----------------------------------------------------------------
	set frame $NotebookPage.corner
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Corner of selection region"
	set cs [$frame childsite]

	foreach plane {x y z} {
	    frame $cs.$plane
	    label $cs.$plane.realvalue \
		    -text "[expr ($Corner($plane)) * $IndexScaling($plane) + $RealBounds(min,$plane)]"

	    scale $cs.$plane.scale \
		    -length 100 \
		    -from $IndexBounds(min,$plane) \
		    -to $IndexBounds(max,$plane) \
		    -variable [scope Corner($plane)] \
		    -orient horizontal -highlightthickness 0

	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $cs.$plane.scale <ButtonRelease> "$this VOIChange"
	}

	#-----------------------------------------------------------------
	# Depth
	#-----------------------------------------------------------------
	set frame $NotebookPage.depth
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Depth of selection region"
	set cs [$frame childsite]

	foreach plane {x y z} {
	    frame $cs.$plane
	    label $cs.$plane.realvalue \
		    -text "[expr $Depth($plane) * $IndexScaling($plane)]"

	    scale $cs.$plane.scale \
		    -length 100 \
		    -from 1 \
		    -to [expr $IndexBounds(max,$plane) - $IndexBounds(min,$plane)+1] \
		    -variable [scope Depth($plane)] \
		    -orient horizontal -highlightthickness 0

	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $cs.$plane.scale <ButtonRelease> "$this VOIChange"
	}


	#-----------------------------------------------------------------
	# Cropping Region
	#-----------------------------------------------------------------
	set frame $NotebookPage.cropping
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Cropping Range"
	set cs [$frame childsite]

	checkbutton $cs.crop \
		-text "Cropping" -variable [scope Crop] \
		-command "$this CropChange"

	foreach plane {min max} {
	    frame $cs.$plane
	    scale $cs.$plane.scale \
		    -length 100 \
		    -from $TableRange(min) \
		    -to $TableRange(max) \
		    -resolution [expr ($TableRange(max) - $TableRange(min)) / 100] \
		    -digits 8 \
		    -variable [scope CroppingRange($plane)] \
		    -orient horizontal -highlightthickness 0 


	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $cs.$plane.scale <ButtonRelease> "$this CroppingRangeChange"
	}

	#-----------------------------------------------------------------
	# Vector Scale Factor
	#-----------------------------------------------------------------
	set frame $NotebookPage.vectorScaleFactor
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Line Length"
	set cs [$frame childsite]

	foreach plane {min max} {
	    frame $cs.$plane
	    scale $cs.$plane.scale \
		    -length 100 \
		    -from 0 \
		    -to 20 \
		    -resolution 0.5 \
		    -variable [scope VectorScaleFactor($plane)] \
		    -orient horizontal -highlightthickness 0 

	    # Need to add button binding here so when slider changes
	    # we update the plane
	    bind $cs.$plane.scale <ButtonRelease> "$this VectorScaleFactorChange"
	}

	#-----------------------------------------------------------------
	# Line Coloring
	#-----------------------------------------------------------------
	set frame $NotebookPage.colorframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Line Color"
	set cs [$frame childsite]

	radiobutton $cs.datacolor \
		-text "Color by data" \
	        -value 1 \
		-variable [scope ColorMap] \
	 	-command "$this ColorMapChange"

	radiobutton $cs.singlecolor \
		-text "Single Color" \
	        -value 0 \
		-variable [scope ColorMap] \
	 	-command "$this ColorMapChange"


	label $cs.vectorlabel -text "Vector Color"
	
	scan $VectorColor "%f %f %f" r g b
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	button $cs.vectorcolor -bg $color \
		-command "$this ExecuteSetVectorColor"

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.visibility
	set cs [$frame childsite]
	pack $cs.visible -anchor w 
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.corner
	set cs [$frame childsite]

	foreach plane {x y z} {
	    pack $cs.$plane.realvalue -side right 
	    pack $cs.$plane.scale -pady 2 -expand yes -fill x 
	    pack $cs.$plane -expand yes -fill x
	}
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.depth
	set cs [$frame childsite]

	foreach plane {x y z} {
	    pack $cs.$plane.realvalue -side right 
	    pack $cs.$plane.scale -pady 2 -expand yes -fill x 
	    pack $cs.$plane -expand yes -fill x
	}
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.cropping
	set cs [$frame childsite]

	pack $cs.crop -anchor w  -pady 4

	foreach plane {min max} {
	    pack $cs.$plane.scale -pady 2 -expand yes -fill x 
	    pack $cs.$plane -expand yes -fill x
	}
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.vectorScaleFactor
	set cs [$frame childsite]

	foreach plane {min max} {
	    pack $cs.$plane.scale -pady 2 -expand yes -fill x 
	    pack $cs.$plane -expand yes -fill x
	}
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.colorframe
	set cs [$frame childsite]
	pack $cs.datacolor -anchor w -pady 4
	pack $cs.singlecolor -anchor w 
	pack $cs.vectorlabel -side left
	pack $cs.vectorcolor -anchor w -pady 4
	pack $cs -fill x
  	pack $frame -fill x

    }   
}

