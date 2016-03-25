##
## File:        cvisAxesInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Axes Interface
##


class cvisAxesInterface {
    inherit cvisInterface

    variable Input ""

    variable AxesType 1

    variable Visible 1

    constructor {name} {
	cvisInterface::constructor $name
    } {
    }

    method Save {file} {
	puts $file "$this SetVisibility $Visible"
	puts $file "$this SetAxesType $AxesType"
    }

    method SetVisibility {visible} {
	if {$Visible != $visible} {
	    set Visible $visible
	    ExecuteSetVisibility
	}
    }
    
    method GetVisibibility {} {
	return $Visible
    }

    method ExecuteSetVisibility {} {
	foreach plane $Input {
	    $plane SetVisibility $Visible
	}

	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetAxesType {type} {
	if {$AxesType != $type} {
	    set AxesType $type
	    ExecuteSetAxesType
	}
    }

    method ExecuteSetAxesType {} {
	foreach plane $Input {
	    $plane SetAxesType $AxesType
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
    }

    method SetNotebookPage {notebook} {
	set NotebookPage $notebook
	#-----------------------------------------------------------------
	# Visibility
	#-----------------------------------------------------------------
	set frame $NotebookPage.visible

	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Visibility"
	set cs [$frame childsite]

	checkbutton $cs.visible \
		-text "Visible" -variable [scope Visible] \
		-command "$this ExecuteSetVisibility"

	#-----------------------------------------------------------------
	# Type
	#-----------------------------------------------------------------
	set frame $NotebookPage.type
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Type"
	set cs [$frame childsite]

	radiobutton $cs.simple -text "Simple" \
		-variable [scope AxesType] \
		-value "1" \
		-command "$this ExecuteSetAxesType"

	radiobutton $cs.float -text "Floating" \
		-variable [scope AxesType] \
		-value "0" \
		-command "$this ExecuteSetAxesType"

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------

	#-----------------------------------------------------------------
	set frame $NotebookPage.visible
	set cs [$frame childsite]
	pack $cs.visible -anchor w 
     	pack $frame -fill x

	#-----------------------------------------------------------------
	set frame $NotebookPage.type
	set cs [$frame childsite]
	pack $cs.simple -side left
	pack $cs.float -side left
  	pack $frame -fill x
    }   
}

