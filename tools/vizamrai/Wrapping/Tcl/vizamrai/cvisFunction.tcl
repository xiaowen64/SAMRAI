##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisFunction.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Slice plane
##

class cvisFunction {
    inherit cvisDisplayObject

    variable Input

    variable Result 0

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
    }

    method SetInput {number input} {
	set Input($number) $input
    }
    
    method GetInput {number} {
	return $Input($number)
    }

    method Function {delta_x delta_y delta_z args} {
	set value [lindex $args 0]
	set Result [expr $Result + $value]
	return 
    }

    method SetResult {value} {
	set Result $value
    }
    
    method GetResult {} {
	return $Result
    }

    method Update {} {
	set array_names [lsort -integer [array names Input]]

	foreach element $array_names {
	    [$Input($element) GetOutput] InitTraversal
	}

	foreach element $array_names {
	    set item($element) [[$Input($element) GetOutput] GetNextItem]
	}
	set base_element $element

	while { [string length $item($base_element)] } {

	    set scaling [$Input($element) GetScaling]

	    scan [$item($base_element) GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)
	    scan [$item($base_element) GetOrigin] "%f %f %f" \
		    origin(x) origin(y) origin(z)
	    foreach plane {x y z} {
		set origin($plane) [expr $origin($plane) / $scaling]
	    }

	    scan [$item($base_element) GetSpacing] "%f %f %f" \
		    delta(x) delta(y) delta(z)
	    foreach plane {x y z} {
		set delta($plane) [expr $delta($plane) / $scaling]
	    }

	    foreach element $array_names {
		set item_scalars($element) \
			[[$item($element) GetPointData] GetScalars]
	    }

	    for {set i 0} {$i < $dim(x) * $dim(y) *$dim(z)} {incr i} {
		# Invoke
		set args ""
		foreach element $array_names {
		    lappend args [$item_scalars($element) GetValue $i]
		}
		eval Function $delta(x) $delta(y) $delta(z) $args
	    }

	    foreach element $array_names {
		set item($element) [[$Input($element) GetOutput] GetNextItem]
	    }
	}
    }
}


