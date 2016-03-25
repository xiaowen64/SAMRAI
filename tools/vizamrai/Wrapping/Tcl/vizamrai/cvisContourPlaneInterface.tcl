##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisContourPlaneInterface.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Slice plane interface
##

class cvisContourPlaneInterface {
    inherit cvisInterface

    common Number 0
    
    variable Input ""

    variable DataIsVisible 0
    
    variable ColorMap 0
    variable ContourColor "1 1 1"
    variable NumberOfContours 10

    variable Jitter 1

    method Save {file} {
	puts $file "$this SetDataVisible $DataIsVisible"
	puts $file "$this SetColorMap $ColorMap"
	puts $file "$this SetContourColor $ContourColor"
	puts $file "$this SetJitter $Jitter"
    }

    constructor {name} {
	cvisInterface::constructor $name
    } {
    }

    method AddInput {sliceplane} {
	set tracervar [::itcl::local cvisTrace #auto]

	lappend Input $sliceplane

	Update

	$sliceplane SetInterface $this
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

    method ColorMapChange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	foreach plane $Input {
	    $plane SetScalarVisibility $ColorMap
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetColorMap {colormap} {
	set tracervar [::itcl::local cvisTrace #auto]
	set ColorMap $colormap
	ColorMapChange
    }

    method SetContourColor {r g b} {
	set ContourColor "$r $g $b"
	
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]

	set frame $NotebookPage.colorframe
	set cs [$frame childsite]
	$cs.contourcolor configure -bg $color
	
	foreach plane $Input {
	    eval $plane SetContourColor $r $g $b
	}
	
	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method ExecuteSetContourColor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	cvisColorChooser color

	eval color SetColor $ContourColor
	
	set value [color Prompt]
	if {[string length $value]} {
	    scan $value "%f %f %f" r g b

	    SetContourColor $r $g $b
	}
	
	delete object color 
    }

    method SetNumberOfContours {number} {
	set tracervar [::itcl::local cvisTrace #auto]
	SetNumberOfContours $NumberOfContours
	ExecuteSetNumberOfContours
    }

    method ExecuteSetNumberOfContours {} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	foreach plane $Input {
	    eval $plane SetNumberOfContours $NumberOfContours
	}

	if {[string length $RenderWindow]} {
	    [$RenderWindow GetRenderWindow] Render
	}
    }

    method SetJitter {jitter} {
	set Jitter $jitter
	ExecuteSetJitter
    }

    method GetJitter {} {
	return $Jitter
    }

    method ExecuteSetJitter {} {
	foreach plane $Input {
	    $plane SetJitter $Jitter
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    # SGS should update interface only, not planes pointed to
    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]
    }

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
	# Contour Lines
	#-----------------------------------------------------------------
	set frame $NotebookPage.contourframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Contour Lines"
	set cs [$frame childsite]

	frame $cs.number
	label $cs.number.numbertitle -text "Number of lines"
	entry $cs.number.number \
		-textvariable [scope NumberOfContours] 
	bind $cs.number.number <Return> "$this ExecuteSetNumberOfContours"

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


	label $cs.contourlabel -text "Contour Color"
	
	scan $ContourColor "%f %f %f" r g b
	set color [format "#%02x%02x%02x" \
		[expr int($r * 255.0)] \
		[expr int($g * 255.0)] \
		[expr int($b * 255.0)]]
	button $cs.contourcolor -bg $color \
		-command "$this ExecuteSetContourColor"

	#-----------------------------------------------------------------
	# Jitter
	#-----------------------------------------------------------------
	set frame $NotebookPage.jitter
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Jitter"
	set cs [$frame childsite]

	scale $cs.scale \
		-from -5 -to 5 \
		-variable [scope Jitter] \
		-orient horizontal -highlightthickness 0

	# Need to add button binding here so when slider changes
	# we update the jitter
	bind $cs.scale <ButtonRelease> "$this ExecuteSetJitter"

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.visibility
	set cs [$frame childsite]
	pack $cs.visible -anchor w 
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.contourframe
	set cs [$frame childsite]
	pack $cs.number.numbertitle -anchor w -side left -pady 8
	pack $cs.number.number -anchor w -side left
	pack $cs.number -anchor w
	pack $cs -fill x
	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.colorframe
	set cs [$frame childsite]
	pack $cs.datacolor -anchor w -pady 4
	pack $cs.singlecolor -anchor w 
	pack $cs.contourlabel -side left
	pack $cs.contourcolor -anchor w -pady 4
	pack $cs -fill x
  	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.jitter
	set cs [$frame childsite]
	pack $cs.scale -anchor w -fill x -pady 4
	pack $cs -fill x
  	pack $frame -fill x
    }   
}


