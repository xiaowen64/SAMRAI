##
## File:        cvisSlicePlaneInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Slice plane interface
##

class cvisSlicePlaneInterface {
    inherit cvisInterface

    common Number 0
    
    variable Slice 0
    variable SlicePlane "z"
    variable Input ""
    variable PlaneInput ""

    variable Scaling 1.0

    variable IndexScaling
    variable IndexBounds

    variable RealBounds

    variable DataIsVisible 1
    variable TextureMap 1
    variable TextureMapInterpolate 0

    variable CellBoundaryWidth 1
    variable CellBoundaryIsVisible 0

    variable NumCellColors 8
    variable CellColors
    variable CellVisible

    variable PatchBoundaryWidth 1
    variable PatchBoundaryIsVisible 0

    variable NumPatchColors 8
    variable PatchColors
    variable PatchVisible

    method Save {file} {
	puts $file "$this SetSlicePlane $SlicePlane"
	puts $file "$this SetSlice $Slice"
	puts $file "$this SetDataVisible $DataIsVisible"
	puts $file "$this SetTextureMap $TextureMap"
	puts $file "$this SetTextureMapInterpolate $TextureMapInterpolate"
	puts $file "$this SetCellBoundaryWidth $CellBoundaryWidth"
	puts $file "$this SetCellBoundaryVisible $CellBoundaryIsVisible"
	set level 0
	while {$level < $NumCellColors} {
	    puts $file "$this SetCellColor $level $CellColors($level)"
	    puts $file "$this SetCellVisibility $level $CellVisible($level)"
	    incr level
	}

	puts $file "$this SetPatchBoundaryWidth $PatchBoundaryWidth"
	puts $file "$this SetPatchBoundaryVisible $PatchBoundaryIsVisible"
	set level 0
	while {$level < $NumPatchColors} {
	    puts $file "$this SetPatchColor $level $PatchColors($level)"
	    puts $file "$this SetPatchVisibility $level $PatchVisible($level)"
	    incr level
	}
    }

    constructor {name} {
	cvisInterface::constructor $name
    } {

	foreach plane {x y z} {
	    set IndexBounds(min,$plane) 0
	    set IndexBounds(max,$plane) 1
	    set IndexScaling($plane) 1.0

	    set RealBounds(min,$plane) 0
	    set RealBounds(max,$plane) 1
	}

	set NumCellColors 8
	# green
	set CellColors(0) "0 1 0"
	# yellow
	set CellColors(1) "1 1 0"
	# light blue
	set CellColors(2) "0 1 1"
	# purple
	set CellColors(3) "1 0 1"
	# pink 
	set CellColors(4) "1 0.75 0.80"
	# orange
	set CellColors(5) "1 0.65 0"
	# red
	set CellColors(6) "1 0 0"
	# blue
	set CellColors(7) "0 0 1"
	
	set CellVisible(0) "1"
	set CellVisible(1) "1"
	set CellVisible(2) "1"
	set CellVisible(3) "1"
	set CellVisible(4) "1"
	set CellVisible(5) "1"
	set CellVisible(6) "1"
	set CellVisible(7) "1"

	set NumPatchColors 8
	# green
	set PatchColors(0) "0 1 0"
	# yellow
	set PatchColors(1) "1 1 0"
	# light blue
	set PatchColors(2) "0 1 1"
	# purple
	set PatchColors(3) "1 0 1"
	# pink 
	set PatchColors(4) "1 0.75 0.80"
	# orange
	set PatchColors(5) "1 0.65 0"
	# red
	set PatchColors(6) "1 0 0"
	# blue
	set PatchColors(7) "0 0 1"
	
	set PatchVisible(0) "1"
	set PatchVisible(1) "1"
	set PatchVisible(2) "1"
	set PatchVisible(3) "1"
	set PatchVisible(4) "1"
	set PatchVisible(5) "1"
	set PatchVisible(6) "1"
	set PatchVisible(7) "1"

    }

    method SetCellColor {level r g b} {
	set CellColors($level) "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]

	set frame $NotebookPage.boundaries
	set cs_main [$frame childsite]
	set frame $cs_main.cellboundary
	set cs [$frame childsite]
	set frame $cs.level_$level
	$frame.level_color configure -bg $color
	
	foreach plane $Input {
	    eval $plane SetCellColor $level $r $g $b
	}
	
	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method ExecuteSetCellColor {level} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $CellColors($level)
	
	set value [color Prompt]
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetCellColor $level $r $g $b
	    
	}
	
	delete object color 
    }

    method SetCellVisibility {level visible} {
	if {$CellVisible($level) != $visible} {
	    set CellVisible($level) $visible
	    ExecuteSetCellVisibility $level
	}
    }

    method GetCellVisibility {level} {
	return $CellVisible($level)
    }

    method ExecuteSetCellVisibility {level} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane $Input {
	    $plane SetCellVisibility $level $CellVisible($level)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetPatchColor {level r g b} {
	set PatchColors($level) "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	
	set frame $NotebookPage.boundaries
	set cs_main [$frame childsite]
	set frame $cs_main.patchboundary
	set cs [$frame childsite]
	set frame $cs.level_$level
	$frame.level_color configure -bg $color

	foreach plane $Input {
	    eval $plane SetPatchColor $level $r $g $b
	}
	
	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method ExecuteSetPatchColor {level} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $PatchColors($level)
	
	set value [color Prompt]
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetPatchColor $level $r $g $b
	    
	}
	
	delete object color 
    }

    method SetPatchVisibility {level visible} {
	if {$PatchVisible($level) != $visible} {
	    set PatchVisible($level) $visible
	    ExecuteSetPatchVisibility $level
	}
    }

    method GetPatchVisibility {level} {
	return $PatchVisible($level)
    }

    method ExecuteSetPatchVisibility {level} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane $Input {
	    $plane SetPatchVisibility $level $PatchVisible($level)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method AddPlaneInput {plane} {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend PlaneInput $plane
    }

    method AddInput {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]

	set Slice 0

	lappend Input $sliceplane

	Update

	$sliceplane SetInterface $this

	SliceChange
    }

    method SetCellBoundaryWidth {width} {
	set tracervar [::itcl::local cvisTrace #auto]
	set CellBoundaryWidth $width
	CellBoundaryWidthChange
    }

    method CellBoundaryWidthChange { } {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetCellBoundaryWidth $CellBoundaryWidth
	}
	
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetCellBoundaryVisible {visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	set CellBoundaryIsVisible $visible
	CellBoundaryVisibleChange
    }

    method CellBoundaryVisibleChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetCellBoundaryVisibility $CellBoundaryIsVisible
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetPatchBoundaryWidth {width} {
	set tracervar [::itcl::local cvisTrace #auto]
	set SetPatchBoundaryWidth $width
	PatchBoundaryWidthChange
    }

    method PatchBoundaryWidthChange { } {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetPatchBoundaryWidth $PatchBoundaryWidth
	}
	
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetPatchBoundaryVisible {visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	set PatchBoundaryIsVisible $visible
	PatchBoundaryVisibleChange
    }

    method PatchBoundaryVisibleChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetPatchBoundaryVisibility $PatchBoundaryIsVisible
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
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

    method TextureMapChange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set frame $NotebookPage.textureframe
	set cs [$frame childsite]
	if {$TextureMap} {
	    $cs.texturemapinterpolate configure -state active
	} {
	    $cs.texturemapinterpolate configure -state disabled
	}
	
	foreach plane $Input {
	    $plane SetTextureMap $TextureMap
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetTextureMap {texturemap} {
	set tracervar [::itcl::local cvisTrace #auto]
	set TextureMap $texturemap
	TextureMapChange
    }

    method TextureMapInterpolateChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    $plane SetTextureMapInterpolate $TextureMapInterpolate
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetTextureMapInterpolate {texturemapinterpolate} {
	set tracervar [::itcl::local cvisTrace #auto]
	set TextureMapInterpolate $texturemapinterpolate
	TextureMapInterpolateChange
    }


    method SetSlice {slice} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Slice $slice

	SliceChange
    }

    method SliceChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	foreach plane $Input {
	    # Put slice in the middle of the cell

	    set real [expr (($Slice + 0.5) * $IndexScaling($SlicePlane)) \
		    + $RealBounds(min,$SlicePlane) ]
	    set screen [expr $real * $Scaling]
	    $plane SetSlice $screen

	    set frame $NotebookPage.scale
	    set cs [$frame childsite]
	    $cs.realvalue configure \
		    -text $real
	}

	foreach plane $PlaneInput {
	    $plane SetSlice $screen
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetSlicePlane {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]
	set SlicePlane $sliceplane
	SlicePlaneChange
    }

    method SlicePlaneChange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	Update

	foreach plane $Input {
	    $plane SetSlicePlane $SlicePlane
	}

	foreach plane $PlaneInput {
	    $plane SetSlicePlane $SlicePlane
	}

	# Plane may need to move
	SliceChange
    }

    method GetSliceInScreenCoord {} {
	set real [expr (($Slice + 0.5) * $IndexScaling($SlicePlane)) \
		+ $RealBounds(min,$SlicePlane) ]
	set screen [expr $real * $Scaling]
	return $screen
    }

    # SGS should update interface only, not planes pointed to
    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]

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

	set real [expr (($Slice + 0.5) * $IndexScaling($SlicePlane)) \
		+ $RealBounds(min,$SlicePlane) ]

	set frame $NotebookPage.scale
	set cs [$frame childsite]

	$cs.scale configure \
		-from $IndexBounds(min,$SlicePlane) -to \
		$IndexBounds(max,$SlicePlane)
	$cs.realvalue configure \
		-text $real
	
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
	# X/Y/Z Plane Selector
	#-----------------------------------------------------------------
	set frame $NotebookPage.plane

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Slice Axis"
	set cs [$frame childsite]


	radiobutton $cs.x -text "X" \
		-variable [scope SlicePlane] \
		-value "x" \
		-command "$this SlicePlaneChange"

	radiobutton $cs.y -text "Y" \
		-variable [scope SlicePlane] \
		-value "y" \
		-command "$this SlicePlaneChange"

	radiobutton $cs.z -text "Z" \
		-variable [scope SlicePlane] \
		-value "z" \
		-command "$this SlicePlaneChange"


	#-----------------------------------------------------------------
	# Slice value
	#-----------------------------------------------------------------
	set frame $NotebookPage.scale
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Slice Location"
	set cs [$frame childsite]

	label $cs.realvalue \
		-text "[expr ($Slice + 0.5) * $IndexScaling($SlicePlane) + $RealBounds(min,$SlicePlane)]"

	scale $cs.scale \
		-length 100 \
		-from $IndexBounds(min,$SlicePlane) \
		-to $IndexBounds(max,$SlicePlane) \
		-variable [scope Slice] \
		-orient horizontal -highlightthickness 0

	# Need to add button binding here so when slider changes
	# we update the plane
	bind $cs.scale <ButtonRelease> "$this SlicePlaneChange"

	#-----------------------------------------------------------------
	# TextureMap 
	#-----------------------------------------------------------------
	set frame $NotebookPage.textureframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Drawing Method"
	set cs [$frame childsite]

	radiobutton $cs.triangles \
		-text "Triangles (slow/accurate)" \
	        -value 0 \
		-variable [scope TextureMap] \
	 	-command "$this TextureMapChange"

	radiobutton $cs.texturemap \
		-text "TextureMap (fast)" \
	        -value 1 \
		-variable [scope TextureMap] \
	 	-command "$this TextureMapChange"

	checkbutton $cs.texturemapinterpolate \
		-text "TextureMap Interpolation" -variable [scope TextureMapInterpolate] \
		-command "$this TextureMapInterpolateChange"

	if {$TextureMap} {
	    $cs.texturemapinterpolate configure -state active
	} {
	    $cs.texturemapinterpolate configure -state disabled
	}

	#-----------------------------------------------------------------
	# Cell Boundary 
	#-----------------------------------------------------------------
	set frame $NotebookPage.boundaries
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Boundary Drawing"
	set cs_main [$frame childsite]

	set frame $cs_main.cellboundary
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Cell"
	set cs [$frame childsite]

	label $cs.edgelabel -text "Edge Width"

	scale $cs.edgewidth \
		-length 20 \
		-from 1 -to 20 \
		-variable [scope CellBoundaryWidth] \
		-orient horizontal -highlightthickness 0

	bind $cs.edgewidth <ButtonRelease> \
		"$this CellBoundaryWidthChange"

	checkbutton $cs.visible \
		-text "Visible" -variable [scope CellBoundaryIsVisible] \
		-command "$this CellBoundaryVisibleChange"

	set level 0
	while {$level < $NumCellColors} {
	    set frame $cs.level_$level
	    frame $frame
	    checkbutton $frame.level_visible \
		    -text "Level $level" -variable [scope CellVisible($level)] \
		    -command "$this ExecuteSetCellVisibility $level"
	    scan $CellColors($level) "%f %f %f" r g b
	    set color [format "#%02x%02x%02x" \
		    [expr int($r * 255.0)] \
		    [expr int($g * 255.0)] \
		    [expr int($b * 255.0)]]
	    button $frame.level_color -bg $color \
		    -command "$this ExecuteSetCellColor $level"
	    incr level
	}

	#-----------------------------------------------------------------
	# Patch Boundary 
	#-----------------------------------------------------------------
	set frame $cs_main.patchboundary
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Patch"
	set cs [$frame childsite]

	label $cs.edgelabel -text "Edge Width"

	scale $cs.edgewidth \
		-length 20 \
		-from 1 -to 20 \
		-variable [scope PatchBoundaryWidth] \
		-orient horizontal -highlightthickness 0

	bind $cs.edgewidth <ButtonRelease> \
		"$this PatchBoundaryWidthChange"

	checkbutton $cs.visible \
		-text "Visible" -variable [scope PatchBoundaryIsVisible] \
		-command "$this PatchBoundaryVisibleChange"

	set level 0
	while {$level < $NumPatchColors} {
	    set frame $cs.level_$level
	    frame $frame
	    checkbutton $frame.level_visible \
		    -text "Level $level" -variable [scope PatchVisible($level)] \
		    -command "$this ExecuteSetPatchVisibility $level"
	    scan $PatchColors($level) "%f %f %f" r g b
	    set color [format "#%02x%02x%02x" \
		    [expr int($r * 255.0)] \
		    [expr int($g * 255.0)] \
		    [expr int($b * 255.0)]]
	    button $frame.level_color -bg $color \
		    -command "$this ExecuteSetPatchColor $level"
	    incr level
	}

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.visibility
	set cs [$frame childsite]
	pack $cs.visible -anchor w 
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.plane
	set cs [$frame childsite]
	pack $cs.x -side left -pady 4
	pack $cs.y -side left
	pack $cs.z -side left
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.scale
	set cs [$frame childsite]
	pack $cs.realvalue -side right 
  	pack $cs.scale -pady 2 -expand yes -fill x 
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.textureframe
	set cs [$frame childsite]
	pack $cs.triangles -anchor w -pady 4
	pack $cs.texturemap -anchor w	
	pack $cs.texturemapinterpolate -anchor e
	pack $cs -fill x
  	pack $frame -fill x

	#-----------------------------------------------------------------

	set frame $NotebookPage.boundaries
	pack $frame -fill x
	set cs_main [$frame childsite]
	pack $cs_main -fill x

	#-----------------------------------------------------------------
	set frame $cs_main.patchboundary
	set cs [$frame childsite]

	pack $frame -pady 4 -side left

	pack $cs.edgelabel
	pack $cs.edgewidth -expand yes -fill x

	pack $cs.visible -anchor w

	set level 0
	while {$level < $NumPatchColors} {
	    set frame $cs.level_$level
	    pack $frame.level_visible -side left
	    pack $frame.level_color -side right
	    incr level
	    pack $frame -anchor w
	}

	#-----------------------------------------------------------------
	set frame $cs_main.cellboundary
	set cs [$frame childsite]

	pack $frame -pady 2 -side left

	pack $cs.edgelabel
	pack $cs.edgewidth -expand yes -fill x

	pack $cs.visible -anchor w

	set level 0
	while {$level < $NumCellColors} {
	    set frame $cs.level_$level
	    pack $frame.level_visible -side left
	    pack $frame.level_color -side right
	    incr level
	    pack $frame -anchor w
	}
    }   
}

