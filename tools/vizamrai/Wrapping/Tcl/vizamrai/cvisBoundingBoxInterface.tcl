##
## File:        cvisBoundingBoxInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Bounding Box Interface
##


class cvisBoundingBoxInterface {
    inherit cvisInterface

    variable Input ""

    variable BoundaryEdgesVisible 1
    variable BoundaryEdgeWidth 1

    variable NumBorderColors 8
    variable BorderColors
    variable BorderVisible

    method Save {file} {
	set level 0
	while {$level < $NumBorderColors} {
	    puts $file "$this SetBorderColor $level $BorderColors($level)"
	    puts $file "$this SetBorderVisibility $level $BorderVisible($level)"
	    incr level
	}
	puts $file "$this SetBoundaryEdgesVisible $BoundaryEdgesVisible"
	puts $file "$this SetBoundaryEdgeWidth $BoundaryEdgeWidth"
    }

    method SelectColor {level} {
	cvisColorChooser color

	eval color SetColor $BorderColors($level)
	
	set value [color Prompt]

	if {[string length $value]} {
	    scan $value "%f %f %f" r g b
	    
	    SetBorderColor $level $r $g $b
	}
	
	delete object color 
    }

    method SetBorderColor {level r g b} {
	set BorderColors($level) "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	
	set frame $NotebookPage.levelcontrol
	set cs [$frame childsite]
	set frame $cs.level_$level
	$frame.level_color configure -bg $color
	
	foreach plane $Input {
	    eval $plane SetBorderColor $level $r $g $b
	}
	
	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method GetBorderColor {$level} {
	return $BorderColors([expr $level % $NumBorderColors])
    }

    constructor {name} {
	cvisInterface::constructor $name
    } {
	set NumBorderColors 8
	# green
	set BorderColors(0) "0 1 0"
	# yellow
	set BorderColors(1) "1 1 0"
	# light blue
	set BorderColors(2) "0 1 1"
	# purple
	set BorderColors(3) "1 0 1"
	# pink 
	set BorderColors(4) "1 0.75 0.80"
	# orange
	set BorderColors(5) "1 0.65 0"
	# red
	set BorderColors(6) "1 0 0"
	# blue
	set BorderColors(7) "0 0 1"
	
	set BorderVisible(0) "1"
	set BorderVisible(1) "1"
	set BorderVisible(2) "1"
	set BorderVisible(3) "1"
	set BorderVisible(4) "1"
	set BorderVisible(5) "1"
	set BorderVisible(6) "1"
	set BorderVisible(7) "1"
    }

    method SetBorderVisibility {level visible} {
	if {$BorderVisible($level) != $visible} {
	    set BorderVisible($level) $visible
	    ExecuteSetBorderVisibility $level
	}
    }

    method GetBorderVisiblity {level} {
	return $BorderVisible($level)
    }

    method ExecuteSetBorderVisibility {level} {
	foreach plane $Input {
	    $plane SetBorderVisibility $level $BorderVisible($level)
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }


    method SetBoundaryEdgesVisible {visible} {
	if {$BoundaryEdgesVisible != $visible} {
	    set BoundaryEdgesVisible $visible
	    ExecuteSetBoundaryEdgesVisible
	}
    }

    method ExecuteSetBoundaryEdgesVisible { } {
	foreach plane $Input {
	    $plane SetVisibility $BoundaryEdgesVisible
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetBoundaryEdgeWidth {width} {
	if {$BoundaryEdgeWidth != $width} {
	    set BoundaryEdgeWidth $width
	    ExecuteSetBoundaryEdgeWidth
	}
    }
    
    method GetBoundaryEdgeWidth {} {
	return $BoundaryEdgeWidth
    }

    method ExecuteSetBoundaryEdgeWidth { } {
	foreach plane $Input {
	    $plane SetBoundaryEdgeWidth $BoundaryEdgeWidth
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method AddInput {sliceplane} {
	lappend Input $sliceplane

	$sliceplane SetInterface $this

	Update
    }

    method Update {} {
	# SGS Sept 5, 2002
	# This was here but I don't see why this is a good thing
	# On an update the width is not changing.  When width 
	# changes this gets called explicitly.

#       ExecuteSetBoundaryEdgeWidth
    }

    method SetNotebookPage {notebook} {
	set NotebookPage $notebook
	#-----------------------------------------------------------------
	# Boundary Edges
	#-----------------------------------------------------------------
	set frame $NotebookPage.visibility

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Visibility"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "Visible" -variable [scope BoundaryEdgesVisible] \
		-command "$this ExecuteSetBoundaryEdgesVisible"

	set frame $NotebookPage.edgewidth

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Edge Width"
	set cs [$frame childsite]

	scale $cs.edgewidth \
		-length 20 \
		-from 1 -to 20 \
		-variable [scope BoundaryEdgeWidth] \
		-orient horizontal -highlightthickness 0

	bind $cs.edgewidth <ButtonRelease> \
		"$this ExecuteSetBoundaryEdgeWidth"

	#-----------------------------------------------------------------
	# Level control
	#-----------------------------------------------------------------

	set frame $NotebookPage.levelcontrol

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Patch Level Visibility and Color"
	set cs [$frame childsite]

	set level 0
	while {$level < $NumBorderColors} {
	    set frame $cs.level_$level
	    frame $frame
	    checkbutton $frame.level_visible \
		    -text "Level $level" -variable [scope BorderVisible($level)] \
		    -command "$this ExecuteSetBorderVisibility $level"
	    scan $BorderColors($level) "%f %f %f" r g b
	    set color [format "#%02x%02x%02x" \
		    [expr int($r * 255.0)] \
		    [expr int($g * 255.0)] \
		    [expr int($b * 255.0)]]
	    button $frame.level_color -bg $color \
		    -command "$this SelectColor $level"
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
	set frame $NotebookPage.edgewidth
	set cs [$frame childsite]
	pack $cs.edgewidth -expand yes -fill x
	pack $frame -fill x

	set frame $NotebookPage.levelcontrol
	pack $frame -fill x
	set cs [$frame childsite]
	set level 0
	while {$level < $NumBorderColors} {
	    set frame $cs.level_$level
	    pack $frame.level_visible -side left
	    pack $frame.level_color -side right
	    incr level
	    pack $frame -anchor w
	}

    }   
}

