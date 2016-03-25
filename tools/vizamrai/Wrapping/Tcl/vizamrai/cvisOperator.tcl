##
## File:        cvisSlicePlane.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Slice plane
##

class cvisOperator {
    inherit cvisDisplayObject

    variable Input

    variable FileName "temp.ascii"

    variable AutoSave 0

    variable TimeProp ""
    variable Time 0.0

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	set TimeProp [cvisProperty $name.TimeProp]
	$name.TimeProp AddCallback "$this UpdateTime"
    }

    method UpdateTime {time} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Time $time
    }

    method GetTime {} {
	set tracervar [::itcl::local cvisTrace #auto]
	set tracervar [::itcl::local cvisTrace #auto]
	return $Time
    }

    method SetTime {time} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Time $time
	$TimeProp SetValue $time
    }

    method GetTimeProp {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$TimeProp GetReference]
    }

    method SetEquivTimeProp { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$TimeProp SetEquiv $property
    }

    method SetInput {number input} {
	set tracervar [::itcl::local cvisTrace #auto]
	set Input($number) $input
    }
    
    method GetInput {number} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Input($number)
    }

    method Function {delta_x delta_y delta_z args} {
	set tracervar [::itcl::local cvisTrace #auto]
	set value [lindex $args 0]
	return $value
    }

    method SetAutoSave { autosave } {
	set tracervar [::itcl::local cvisTrace #auto]
	set AutoSave $autosave
    }

    method GetAutoSave { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $AutoSave
    }

    method GetFileName {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $FileName
    }

    method SetFileName {filename} {
	set tracervar [::itcl::local cvisTrace #auto]
	set FileName $filename
    }

    method SaveFileAsASCII {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set file [open "$FileName" w]

	puts $file $Time

	foreach element [array names Input] {
	    [$Input($element) GetOutput] InitTraversal
	}

	# Write out the number of patches
	puts $file [[$Input($element) GetOutput] GetNumberOfItems]

	foreach element [array names Input] {
	    set item($element) [[$Input($element) GetOutput] GetNextItem]
	}
	set base_element $element

	while { [string length $item($base_element)] } {

	    set scaling [$Input($element) GetScaling]

	    # Write out the number dimensions of this patch
	    scan [$item($base_element) GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)
	    puts $file "$dim(x) $dim(y) $dim(z)"

	    # Write out the origin for this patch
	    scan [$item($base_element) GetOrigin] "%f %f %f" \
		    origin(x) origin(y) origin(z)
	    foreach plane {x y z} {
		set origin($plane) [expr $origin($plane) / $scaling]
	    }
	    puts $file "$origin(x) $origin(y) $origin(z)"

	    # Write out the delta spacing for this patch
	    scan [$item($base_element) GetSpacing] "%f %f %f" \
		    delta(x) delta(y) delta(z)
	    foreach plane {x y z} {
		set delta($plane) [expr $delta($plane) / $scaling]
	    }
	    puts $file "$delta(x) $delta(y) $delta(z)"

	    foreach element [array names Input] {
		set item_scalars($element) \
			[[$item($element) GetPointData] GetArray 0]
	    }

	    for {set i 0} {$i < $dim(x) * $dim(y) *$dim(z)} {incr i} {
		# Invoke
		set args ""
		foreach element [array names Input] {
		    lappend args [$item_scalars($element) GetValue $i]
		}
		set value [eval Function $delta(x) $delta(y) $delta(z) $args]
		puts $file $value
	    }

	    foreach element [array names Input] {
		set item($element) [[$Input($element) GetOutput] GetNextItem]
	    }
	}

	close $file
    }

    method Update {} {
	if {$AutoSave} {
	    ::cvis::queue::Add "$this SaveFileAsASCII"
	}
    }
}


