##
## File:        cvisIsoInterface.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Isosurface Interface
##

class cvisIsoInterface {
    inherit cvisInterface

    common Number 0
    
    variable IsoValue 0.5
    variable TableRange 

    variable Input ""

    constructor {name} {
	cvisInterface::constructor $name
    } {
	set TableRange(min) 0.0
	set TableRange(max) 1.0
    }

    method Save {file} {
	puts $file "$this SetIsoValue $IsoValue"
    }

    method AddInput {iso} {
	lappend Input $iso

	$iso SetInterface $this
    }

    method ExecuteSetIsoValue {} {
	foreach item $Input {
	    $item SetIsoValue $IsoValue
	}
	if {[string length $RenderWindow]} {
	    $RenderWindow Render
	}
    }

    method SetIsoValue {isovalue} {
	if {$IsoValue != $isovalue} {
	    set IsoValue $isovalue
	    ExecuteSetIsoValue
	}
    }

    method SetTableRange {min max} {
	set TableRange(min) $min
	set TableRange(max) $max

	ExecuteSetTableRange
    }

    method GetTableRange {} {
	return "$TableRange(min) $TableRange(max)"
    }

    method ExecuteSetTableRange {} {
	if { $IsoValue < $TableRange(min) } {
	    set IsoValue $TableRange(min)
	}

	if { $IsoValue > $TableRange(max) } {
	    set IsoValue $TableRange(max)
	}

	set frame $NotebookPage.isovalueframe
	set cs [$frame childsite]

	$cs.scale configure \
		-from $TableRange(min) -to $TableRange(max) \
		-resolution [expr ($TableRange(max) - $TableRange(min)) / 100] \
		-digits 8
    }

    method Update {} {
    }

    method SetNotebookPage {notebook} {
	set NotebookPage $notebook

	#-----------------------------------------------------------------
	# IsoValue value
	#-----------------------------------------------------------------

	set frame $NotebookPage.isovalueframe
	iwidgets::Labeledframe $frame \
		-labelpos $::cvis::gui::PanelHeaderPosition \
		-labeltext "Iso Surface"
	set cs [$frame childsite]

	scale $cs.scale \
		-length 100 \
		-from $TableRange(min) -to $TableRange(max) \
		-resolution [expr ($TableRange(max) - $TableRange(min)) / 100] \
		-digits 3 \
		-variable [scope IsoValue] \
		-orient horizontal -highlightthickness 0

	# Need to add button binding here so when slider changes
	# we update the plane
	bind $cs.scale <ButtonRelease> \
		"$this ExecuteSetIsoValue"

	#-----------------------------------------------------------------
	# Pack widgets
	#-----------------------------------------------------------------
	set frame $NotebookPage.isovalueframe
	set cs [$frame childsite]

  	pack $cs.scale -side left -expand yes -fill x -pady 6
	
  	pack $cs -fill x
  	pack $frame -fill x
    }   
}

