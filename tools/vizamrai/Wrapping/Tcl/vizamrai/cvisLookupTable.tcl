##
## File:        cvisLookupTable.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Color Lookup Table
##

class cvisLookupTable { 
    inherit cvisDisplayObject

    variable LookupTable ""

    variable Scale 0
    variable Ramp 0

    variable TableRange "0 1"
    variable AutoTableRangeScaling 1
    
    variable HueRange "0.66667 0.0"

    variable Visible 1

    variable Reader ""

    variable K 5.0
    variable K0 -0.001
    variable K1 0.5

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {




	$this.colorBar SetLookupTable $this
    }

    method SetRenderWindow {renderwindow} {
	set RenderWindow $renderwindow
	$this.colorBar SetRenderWindow $renderwindow
    }

    variable ColorFileName
    method SetColorFileName { filename } {
	set ColorFileName $filename
	if {[string length $filename]} { 
	    set fp [open $filename r]
	    if { $fp == "" } {
		puts "Error opening colormap file <$filename>"
	    } else {
		for {set index 0} {$index < 256} {incr index} {
		    gets $fp line
		    scan $line "%i %i %i" r g b
			    [expr $r/256.0] \
			    [expr $g/256.0] \
			    [expr $b/256.0] \
			    1.0
		}
	    }
	    # Update things that need to know about changes to 
	    # colormaps
	    UpdateTableRange
	}
    }

    method GetColorFileName {} {
	return $ColorFileName
    }

    variable AlphaFileName
    method SetAlphaFileName { filename } {
	set AlphaFileName $filename
	if {[string length $filename]} { 
	    set fp [open $filename r]
	    if { $fp == "" } {
		puts "Error opening alpha map file <$filename>"
	    } else {
		for {set index 0} {$index < 256} {incr index} {
		    gets $fp line
		    scan $line "%i %i %i" alpha a2 a3
			    [lindex $oldvalues 0] \
			    [lindex $oldvalues 1] \
			    [lindex $oldvalues 2] \
			    [expr $alpha/256.0]
		}
	    }
	    # Update things that need to know about changes to 
	    # colormaps
	    UpdateTableRange
	}
    }

    method GetAlphaFileName {} {
	return $AlphaFileName
    }

    method SetType {type} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$Scale != $type} {
	    set Scale $type
	    if {[string length $LookupTable]} {
		$LookupTable SetScale $Scale
#		eval $LookupTable SetTableRange $TableRange
		eval $LookupTable SetHueRange $HueRange

		UpdateTableRange

		foreach item $LogDepends {
		    eval $item SetType $Type
		}
	    }
	}
	UpdateLookupTable
    }

    method GetType {} {
	return $Type
    }

    method SetK {k} {
	if {$k != $K} {
	    set K $k
	    $LookupTable SetK $K
	    UpdateLookupTable
	}
    }
    method GetK {} {
	return $K
    }

    
    method GetK0 {} {
	return $K0
    }
    method SetK0 {k} {
	if {$k != $K0} {
	    set K0 $k
	    $LookupTable SetK0 $K0
	    UpdateLookupTable
	}
    }
    
    method GetK1 {} {
	return $K1
    }
    method SetK1 {k} {
	if {$k != $K1} {
	    set K1 $k
	    $LookupTable SetK1 $K1
	    UpdateLookupTable
	}
    }

    method SetRamp {ramp} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$Ramp != $ramp} {
	    set Ramp $ramp
	    if {[string length $LookupTable]} {
		$LookupTable SetRamp $Ramp
		eval $LookupTable SetHueRange $HueRange


		UpdateTableRange
	    }
	}
	UpdateLookupTable
    }

    method GetRamp {} {
	return $Ramp
    }

    method SetVisibility {visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$Visible != $visible} {
	    set Visible $visible
	    $this.colorBar SetVisibility $visible
	}
    }

    method SetTableRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]
	set TableRange "$min $max"
	if {!$AutoTableRangeScaling} {
	    UpdateTableRange
	}
    }

    method SetAutoTableRangeScaling {autoscale} {
	set tracervar [::itcl::local cvisTrace #auto]
	if { $AutoTableRangeScaling != $autoscale } {
	    set AutoTableRangeScaling $autoscale
	    UpdateTableRange
	}
    }

    method UpdateTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if { $AutoTableRangeScaling } {
	    eval $LookupTable SetTableRange [$Reader GetScalarRange]
	    foreach item $Depends {
		eval $item SetTableRange [$Reader GetScalarRange]
	    }
	} {
	    eval $LookupTable SetTableRange $TableRange
	    foreach item $Depends {
		eval $item SetTableRange $TableRange
	    }
	}
    }

    method SetInput { reader } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Reader $reader
	$reader AddOutput $this
    }
    
    method GetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $TableRange
    }

    method SetHueRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]
	set HueRange "$min $max"
	$LookupTable SetHueRange $min $max
	UpdateLookupTable
    }

    method GetColorBar {} {
	return $this.colorBar
    }

    method GetLookupTable { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $LookupTable
    }

    method UpdateLookupTable {} {
	set tracervar [::itcl::local cvisTrace #auto]

	$this.colorBar Update

	# Update things that depend on this
	foreach item $Depends {
	    $item UpdateLookupTable
	}
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$AutoTableRangeScaling} {
	    UpdateTableRange
	}
    }

    variable Depends ""
    method AddDepends { depend } {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Depends $depend
    }

    variable LogDepends ""
    method AddLogDepends { depend } {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend LogDepends $depend
    }

    method AddOutput {list} {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Outputs $list
    }
}
