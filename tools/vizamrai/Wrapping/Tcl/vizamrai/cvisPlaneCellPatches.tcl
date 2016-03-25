##
## File:        cvisPlaneCellPatches
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Build a set of patches which extend cell data to the
##              patch boundaries for use in carpet plots
##

class cvisPlaneCellPatches {
    inherit cvisDisplayObject

    variable Reader ""

    variable SlicePlane "z"
    variable Slice 0.5

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
    }

    method SetScalarRange {min max} {
	set ScalarRange "$min $max"
    }

    method GetScalarRange {} {
	return $ScalarRange 
    }

    method SetInput { reader } {
	set Reader $reader
	$reader AddOutput $this
    }

    method GetBounds { } { 
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetBounds]
	} {
	    return {0.0 1.0 0.0 1.0 0.0 1.0}
	}
    }

    method GetSpacing { } {
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetSpacing] 
	} {
	    return {0.01 0.01 0.01}
	}
    }

    method SetSlice {value} {
	if {[expr $value != $Slice]} {
	    set Slice $value
	    Execute
	}
    }

    method GetSlice {} {
	return $Slice
    }

    method GetIntersectingBoxes {} {

	set intersectedboxes ""

	set collection [[$Reader GetReader] GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {

	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    real(min,x) real(max,x) \
		    real(min,y) real(max,y) \
		    real(min,z) real(max,z)

	    puts "$this Actually Slicing $Slice"

	    if { $Slice >= $real(min,$SlicePlane) && \
		    $Slice <= $real(max,$SlicePlane) } {
		lappend intersectedboxes $patch
	    }
	    
	    set item [$collection GetNextItem]
	    incr patch
	}

	return $intersectedboxes
    }

    method SortedIntersectionList {plane patchlist lowerintersectarray upperintersectarray } {
	upvar $lowerintersectarray lowerlist
	upvar $upperintersectarray upperlist

	switch $SlicePlane {
	    "z" {
		switch $plane {
		    "x" {
			set intersect_axis y
			set proj_axis x
		    }
		    "y" {
			set intersect_axis x
			set proj_axis y
		    }
		}
	    }
	}
	
	foreach patch $patchlist {
	    
	    set box [$Reader GetPatch $patch]

	    scan [$box GetBounds] "%f %f %f %f %f %f" \
		    real(min,x) real(max,x) \
		    real(min,y) real(max,y) \
		    real(min,z) real(max,z)

	    lappend t_lowerlist($real(min,$proj_axis)) \
		    "$real(min,$intersect_axis) $patch $real(max,$intersect_axis) $patch"

	    lappend t_upperlist($real(max,$proj_axis)) \
		    "$real(min,$intersect_axis) $patch $real(max,$intersect_axis) $patch"

	    set intersectlines($real(min,$proj_axis)) 1
	    set intersectlines($real(max,$proj_axis)) 1
	}

	foreach i [array names t_lowerlist] {
	    set t_lowerlist($i) [lsort -real -index 0 $t_lowerlist($i)]
	    foreach element $t_lowerlist($i) {
		lappend lowerlist($i) "[lindex $element 0] [lindex $element 1]"
		lappend lowerlist($i) "[lindex $element 2] [lindex $element 3]"
	    }
	}

	foreach i [array names t_upperlist] {
	    set t_upperlist($i) [lsort -real -index 0 $t_upperlist($i)]
	    foreach element $t_upperlist($i) {
		lappend upperlist($i) "[lindex $element 0] [lindex $element 1]"
		lappend upperlist($i) "[lindex $element 2] [lindex $element 3]"
	    }
	}

	return [lsort [array names intersectlines]]
    }

    method InsertEmptyBoxes { _intersectarray } {
	upvar $_intersectarray intersectarray

	foreach i [array names intersectarray] {
	    
	    set list $intersectarray($i)
	    set new_list ""

	    set top_point [lindex $list 0]
	    while {[llength $list]} {
		set bottom_point [lindex $list 0]
		set list [lreplace $list 0 0]

		if { [lindex $top_point 0] == [lindex $bottom_point 0] } {
		    lappend new_list $top_point
		    lappend new_list $bottom_point
		} {
		    lappend new_list $top_point
		    lappend new_list "[lindex $top_point 0] -1"
		    lappend new_list "[lindex $bottom_point 0] -1"
		    lappend new_list $bottom_point
		}

		set top_point [lindex $list 0]
		set list [lreplace $list 0 0]
	    }
	    
	    # strip of duplicate of first point (was the starting bottom)
	    set new_list [lreplace $new_list 0 0]
	    lappend new_list "[lindex $top_point 0] -1"

	    set intersectarray($i) $new_list
	}


    }

    variable NumberOfPatches 0

    method CreateTriangleStrips {intersects listarray} {
	upvar $listarray list

	set data [$Reader GetReader]


	    incr NumberOfPatches -1
	    while {$NumberOfPatches > -1} { 
		incr NumberOfPatches -1
	    }
	    set NumberOfPatches 0
	}


	# SGS how do we compute this?
	set point 0



	
#	foreach i [lindex $intersects 0] 
	foreach i $intersects {
	    while {[llength $list($i)]} {
		set line_point 0

		set start_point [lindex $list($i) 0]
		set list($i) [lreplace $list($i) 0 0]

		set end_point [lindex $list($i) 0]
		set list($i) [lreplace $list($i) 0 0]

		set starting_value [lindex $start_point 0] 
		set ending_value [lindex $end_point 0] 

		set left_patch [lindex $start_point 1]
		set right_patch [lindex $start_point 2]
		
		# If not an empty space
		if {$left_patch >= 0} {
		    set left_level [$data GetPatchLevel $left_patch]
		    set left_data [$data GetPatch $left_patch]
		}

		# If not an empty space
		if {$right_patch >= 0} {
		    set right_level [$data GetPatchLevel $right_patch]
		    set right_data [$data GetPatch $right_patch]
		}

		
		if {$left_patch < 0} {

		    set right_level [$data GetPatchLevel $right_patch]

		    scan [$right_data GetSpacing] "%f %f %f" \
			    spacing(x) spacing(y)  spacing(z)

		    set dim(y) [expr  \
			    round(($ending_value - $starting_value) / $spacing(y))]


		    
			    SetNumberOfIds [expr $dim(y) * 2 + 2]

		    #Insert 0'th primer point
			    $i \
			    $starting_value \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point

			    $i \
			    [expr $starting_value + 0.5 * $spacing(y)] \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point


			    [expr $i + 0.5 * $spacing(x)] \
			    [expr $starting_value + 0.5 * $spacing(y)] \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point


			    $i \
			    [expr $starting_value + 1.0 * $spacing(y)] \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point

		    for {set index 1} {$index < $dim(y)} \
			    {incr index} {
				[expr $i + 0.5 * $spacing(x)] \
				[expr $starting_value+($index+0.5)*$spacing(y)]\
				$Slice
				SetId $line_point $point

			incr point
			incr line_point

				$i \
				[expr $starting_value+($index+1.0)*$spacing(y)]\
				$Slice
				SetId $line_point $point

			incr point
			incr line_point

		    }

		    
		    incr NumberOfPatches

		} elseif {$right_patch < 0} {
		    set left_level [$data GetPatchLevel $left_patch]

		    scan [$left_data GetSpacing] "%f %f %f" \
			    spacing(x) spacing(y)  spacing(z)

		    set dim(y) [expr  \
			    round(($ending_value - $starting_value) / $spacing(y))]

		    
			    SetNumberOfIds [expr $dim(y) * 2 + 2]
		    
		    #Insert 0'th primer point
			    $i \
			    $starting_value \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point

			    $i \
			    [expr $starting_value + 0.5 * $spacing(y)] \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point


			    [expr $i - 0.5 * $spacing(x)] \
			    [expr $starting_value + 0.5 * $spacing(y)] \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point


			    $i \
			    [expr $starting_value + 1.0 * $spacing(y)] \
			    $Slice
			    SetId $line_point $point
		    incr point
		    incr line_point

		    for {set index 1} {$index < $dim(y)} \
			    {incr index} {
				[expr $i - 0.5 * $spacing(x)] \
				[expr $starting_value+($index+0.5)*$spacing(y)]\
				$Slice
				SetId $line_point $point

			incr point
			incr line_point

				$i \
				[expr $starting_value+($index+1.0)*$spacing(y)]\
				$Slice
				SetId $line_point $point

			incr point
			incr line_point
		    }

		    
		    incr NumberOfPatches
		} else {
		    # Need to do both
		}
	    }

	}




    }

    method StripList {intersects lowerintersectarray upperintersectarray listarray} {
	upvar $lowerintersectarray lowerlist
	upvar $upperintersectarray upperlist
	upvar $listarray list

	set boxes ""

	foreach i $intersects {

	    if {[info exists lowerlist($i)]} {
		set lower $lowerlist($i)
	    } {
		set lower ""
	    }
	    lappend lower "$::cvis::math::FLT_MAX -1"
	    set currentlowerbox -1

	    if {[info exists upperlist($i)]} {
		set upper $upperlist($i)

	    } {
		set upper ""
	    }
	    lappend upper "$::cvis::math::FLT_MAX -1"
	    set currentupperbox -1

	    set currentlocation $::cvis::math::FLT_MIN

	    set list($i) ""
	    lappend list($i) "$currentlocation $currentupperbox $currentlowerbox"

	    while { $currentlocation < $::cvis::math::FLT_MAX } {
		set lower_temp [lindex $lower 0]
		set upper_temp [lindex $upper 0]
		set lower_location [lindex $lower_temp 0]
		set upper_location [lindex $upper_temp 0]

		if { $lower_location > $upper_location } {
		    set currentlocation $upper_location
		    set currentupperbox [lindex $upper_temp 1]
		    set upper [lreplace $upper 0 0]
		} elseif { $lower_location < $upper_location} {
		    set currentlocation $lower_location
		    set currentlowerbox [lindex $lower_temp 1]
		    set lower [lreplace $lower 0 0]
		} {
		    set currentlocation $lower_location

		    set currentupperbox [lindex $upper_temp 1]
		    if { $upper_location < $::cvis::math::FLT_MAX } {

			set upper [lreplace $upper 0 0]
		    }

		    set currentlowerbox [lindex $lower_temp 1]
		    if { $lower_location < $::cvis::math::FLT_MAX } {
			set lower [lreplace $lower 0 0]
		    }
		}

		lappend list($i) "$currentlocation $currentupperbox $currentlowerbox"
	    }

	    set list($i) [lreplace $list($i) 0 0]
	    set list($i) [lreplace $list($i) end end]
	}
    }

    method Execute {} {
	if {![string length $RenderWindow] } { 
	    return
	}

	# First Get a list of boxes that intersect with the current
	# slice plane
	set intersectedboxes [GetIntersectingBoxes]

	set intersects [SortedIntersectionList "x" $intersectedboxes \
		leftintersections rightintersections]

	InsertEmptyBoxes leftintersections 
	InsertEmptyBoxes rightintersections

	StripList $intersects leftintersections rightintersections striplist

	CreateTriangleStrips $intersects striplist

    }

    method SetSlicePlane {plane} {
	set SlicePlane $plane

	set Create 1
	Update
    }
    variable Depends ""

    method AddDepends { depend } {
	lappend Depends $depend
    }

    method Update {} {

	if $Create { 
	}

	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}
    }
}

