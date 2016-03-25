##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisContourPlane.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Slice plane
##

class cvisContourPlane {
    inherit cvisDisplayObject

    # Slice plane
    variable SlicePlane "z"
    variable Slice 0

    # Data slice
    variable DataExists 0
    variable DataIsVisible 0

    # Contouring 
    # controls coloring of contour lines
    variable ScalarVisibility 0
    variable ContourColor "1 1 1"
    variable NumberOfContours 10
    variable Jitter 1


    variable TableRange "0.0 1.0"
    variable AxisScale "1.0 1.0 1.0"

    # Links to other objects
    variable Reader ""
    variable Interface ""
    variable LookupTable ""
    variable Depends ""
    variable PreviousCollection ""

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
    }

    #=======================================================================
    # Access Functions
    #=======================================================================

    method GetScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [$Reader GetScaling]
	} {
	    return 1.0
	}
    }

    method GetBounds { } { 
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetBounds]
	} {
	    return {0.0 1.0 0.0 1.0 0.0 1.0}
	}
    }


    method GetIndexScaling { } { 
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetIndexScaling]
	} {
	    return {1.0 1.0 1.0}
	}
    }

    method GetIndexBounds { } { 
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetIndexBounds]
	} {
	    return {0 1 0 1 0 1}
	}
    }

    method GetSpacing { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetSpacing] 
	} {
	    return {0.01 0.01 0.01}
	}
    }

    method SetSlice {value} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$Slice == $value} {
	    return
	}

	set Slice $value
	::cvis::queue::Add "$this ExecuteSetSlice"
	
    }

    method ExecuteSetSlice {} {
	ExecuteDataSetSlice
    }

    method SetSlicePlane {plane} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$SlicePlane == $plane} {
	    return
	}

	set SlicePlane $plane

	Update
    }

    method SetContourColor { r g b} {
	set tracervar [::itcl::local cvisTrace #auto]

	set ContourColor "$r $g $b"
	
	if {[string length $Reader]} {
	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		if { [string length [info commands \
			$this.$item.planeActor]] } {
		    eval [$this.$item.planeActor GetProperty] SetColor $ContourColor
		}
		set item [$collection GetNextItem]
	    }
	}
    }

    method GetContourColor {} {
	return $ContourColor
    }

    method SetNumberOfContours { number } {
	if {$NumberOfContours != $number} {
	    set NumberOfContours $number
	    if {[string length $Reader]} {
		set collection [$Reader GetCollection]
		$collection InitTraversal
		set item [$collection GetNextItem]
		while {[string length $item]} {
		    if { [string length [info commands \
		    }
		    set item [$collection GetNextItem]
		}
	    }
	}
    }

    method GetNumberOfContours {} {
	return $NumberOfContours
    }

    method SetScalarVisibility {mode} {
	if {$ScalarVisibility != $mode } {
	    set ScalarVisibility $mode

	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		set item [$collection GetNextItem]
	    }
	}
    }

    method GetScalarVisibilty {} {
	return $ScalarVisibility
    }

    method SetJitter { jitter } {
	if {$Jitter == $jitter} {
	    return
	}

	set Jitter $jitter

	if {$DataExists} {
	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]

	    while {[string length $item]} {
		lappend PreviousCollection $item

		foreach plane {x y z} i {0 1 2}  {
		    set origin($plane)  [lindex [$item GetOrigin] $i]
		    set spacing($plane) [lindex [$item GetSpacing] $i]   
		}
		
		foreach plane {x y z} {
		    set origin_new($plane) [expr $origin($plane) + \
			    $spacing($plane)/2.0]

		}
		
		set origin_new($SlicePlane) [expr $origin_new($SlicePlane) + \
			double($Jitter) * $spacing($SlicePlane) * \
			$::cvis::math::JITTER]
		
			$origin_new(x) \
			$origin_new(y) \
			$origin_new(z)

		set item [$collection GetNextItem]
	    }
	}
    }

    method GetJitter {} {
	return $Jitter
    }

    #=======================================================================
    # Data Slice
    #=======================================================================
    method SetDataVisibility {vis} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$vis == $DataIsVisible} {
	    return
	}
	
	if {$DataIsVisible} {
	    set DataIsVisible 0
	    ::cvis::queue::Add "$this ExecuteDataDelete"
	} {
	    set DataIsVisible 1
	    ::cvis::queue::Add "$this ExecuteDataCreate"
	}
    }
    
    method GetDataVisibility {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $DataIsVisible
    }

    method ExecuteDataCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	if {$DataExists} {
	    return
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}
	
	set DataExists 1
	
	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]

	while {[string length $item]} {
	    lappend PreviousCollection $item

	    foreach plane {x y z} i {0 1 2}  {
		set origin($plane) [lindex [$item GetOrigin] $i]
		set dim($plane)    [lindex [$item GetDimensions] $i]
		set spacing($plane)    [lindex [$item GetSpacing] $i]
	    }

	    foreach i {0 1 2 3 4 5 } {
		set bounds($i) [lindex [$item GetBounds] $i]
	    }

	    set IndexSlice [expr int(($Slice-$origin($SlicePlane)) / \
		    $spacing($SlicePlane))]

	    foreach plane "x y z" {
		if { $dim($plane) > 1 } {
		    set newdata_dim($plane) [expr $dim($plane) - 1]
		} {
		    set newdata_dim($plane) $dim($plane)
		}

	    }

	    # Convert to Cell to Node centered data
		    $newdata_dim(x) $newdata_dim(y) $newdata_dim(z)

		    [[$item GetCellData] GetScalars]

	    foreach plane {x y z} {
		set origin_new($plane) [expr $origin($plane) + \
			$spacing($plane)/2.0]
		
	    }

	    set origin_new($SlicePlane) [expr $origin_new($SlicePlane) + \
		    double($Jitter) * $spacing($SlicePlane) * \
		    $::cvis::math::JITTER]

		    $origin_new(x) \
		    $origin_new(y) \
		    $origin_new(z)
		    
		    $spacing(x) \
		    $spacing(y) \
		    $spacing(z) 



	    switch $SlicePlane {
		"x" {
			    $IndexSlice $IndexSlice \
			    0 $dim(y) \
			    0 $dim(z) 
		}
		"y" {
			    0 $dim(x) \
			    $IndexSlice $IndexSlice \
			    0 $dim(z) 
		}
		"z" {
			    0 $dim(x) \
			    0 $dim(y) \
			    $IndexSlice $IndexSlice
		}
	    }


	    # SGS need to generate this from the GUI

	    
	    if {[string length $LookupTable] } { 
			[$LookupTable GetLookupTable]
	    }
	    
	    cvisActor $this.$item.planeActor

	    
	    scan [$this.$item.planeActor GetPosition] "%f %f %f" \
		    position(x) position(y) position(z)

	    $this.$item.planeActor SetPosition \
		    $position(x) $position(y) $position(z)

	    # Turn Actor off by default since it might not lie on the
	    # current slice plane.
	    $this.$item.planeActor SetVisibility 0
	    
	    $RenderWindow AddActor $this.$item.planeActor

	    set item [$collection GetNextItem]
	}
	
	$RenderWindow SetScaling
	
	ExecuteDataSetSlice 
    }

    method ExecuteDataDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set DataExists 0
	foreach item $PreviousCollection {

	    if {[string length \
	    }

	    if {[string length \
		$RenderWindow RemoveActor $this.$item.planeActor
		::itcl::delete object $this.$item.planeActor
	    }
	}
    }

    # SGS this method is broken?
    method ExecuteDataSetSlice {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}

	# If no data reader then can't do anything
	if {![string length $Reader]} {
	    return
	}

	# If we have not gotten an Update from reader then there
	# is nothing to do.
	if {!$DataExists} {
	    return
	}

	scan [[$Reader GetReader] GetBounds] "%f %f %f %f %f %f" \
		real(min,x) real(max,x) \
		real(min,y) real(max,y) \
		real(min,z) real(max,z)

	# Push the slice in from the boundaries
	set jitter [expr ($real(max,$SlicePlane) - $real(min,$SlicePlane))\
		* $::cvis::math::BOUNDARY_JITTER ]
	
	if { $Slice < [expr $real(min,$SlicePlane) + $jitter] } {
	    set Slice [expr $Slice + $jitter]
	}
	
	if { $Slice > [expr $real(max,$SlicePlane) - $jitter] } {
	    set Slice [expr $Slice - $jitter]
	}
	
	set collection [[$Reader GetReader] GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    
	    # Need to do more complicated stuff here for 
	    # see the vizamrai for more information
	    # need to set to ends and set vizibility off
	    
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    real(min,x) real(max,x) \
		    real(min,y) real(max,y) \
		    real(min,z) real(max,z)
	    
	    if { $Slice >= $real(min,$SlicePlane) && \
		    $Slice <= $real(max,$SlicePlane) } {

		foreach plane {x y z} i {0 1 2}  {
		    set origin($plane) [lindex [$item GetOrigin] $i]
		    set dim($plane)    [lindex [$item GetDimensions] $i]
		    set spacing($plane)    [lindex [$item GetSpacing] $i]
		}
		
		foreach i {0 1 2 3 4 5 } {
		    set bounds($i) [lindex [$item GetBounds] $i]
		}
		
		set IndexSlice [expr int(($Slice-$origin($SlicePlane)) / \
			$spacing($SlicePlane))]

		switch $SlicePlane {
		    "x" {
				$IndexSlice $IndexSlice \
				0 $dim(y) \
				0 $dim(z) 
		    }
		    "y" {
				0 $dim(x) \
				$IndexSlice $IndexSlice \
				0 $dim(z) 
		    }
		    "z" {
				0 $dim(x) \
				0 $dim(y) \
				$IndexSlice $IndexSlice
			
			if {$dim(x) == 2} {
			}
		    }
		}
		
		$this.$item.planeActor SetVisibility $DataIsVisible
		$this.$item.planeActor Update
	    } {
		$this.$item.planeActor SetVisibility 0
	    }
	    
	    set item [$collection GetNextItem]
	}
    }

    method SetAxisScale { scalex scaley scalez } {
	set tracervar [::itcl::local cvisTrace #auto]
	set AxisScale "$scalex $scaley $scalez"
	
	set collection [[$Reader GetReader] GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    eval $this.$item.planeActor SetScale $AxisScale
	    eval $this.$item.iedgeActor SetScale $AxisScale
	    set item [$collection GetNextItem]
	}
    }

    method GetAxisScale { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $AxisScale
    }

    method SetInput { reader } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Reader $reader
	$reader AddOutput $this
    }

    method SetInterface { interface } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Interface $interface
	$interface Update
    }

    method GetSlice {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Slice
    }

    method SetRenderWindow {window} {
	set tracervar [::itcl::local cvisTrace #auto]
	set RenderWindow $window
    }
    
    method SetTableRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]

	set TableRange "$min $max"

	::cvis::queue::Add "$this ExecuteSetTableRange"
    }

    method ExecuteSetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	    while {[string length $item]} {
		set item [$collection GetNextItem]
	    }
	}
    }

    method GetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $TableRange 
    }

    method SetLookupTable {lut} {
	set tracervar [::itcl::local cvisTrace #auto]
	$lut AddDepends $this
	set LookupTable $lut
	UpdateLookupTable
    }
    
    method UpdateLookupTable {} {
	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    if {[string length $LookupTable] } { 
			    [$LookupTable GetLookupTable]
		}
	    }
	    set item [$collection GetNextItem]
	}

    }

    method AddDepends { depend } {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Depends $depend
    }

    method ExecuteSetPreviousCollection {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set PreviousCollection ""
	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    lappend PreviousCollection $item
	    set item [$collection GetNextItem]
	}
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]

	# Update user interface elements
	if {[string length $Interface]} {
	    $Interface Update
	    # Get values that might have changed based on UI control
	    # SGS this is bad don't have this anymore
	    # set Slice [$Interface GetSliceInScreenCoord]
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	# Schedule CPU intensive tasks that accomplish the update
	if {$DataExists} {
	    set DataExists 0
	    ::cvis::queue::Add "$this ExecuteDataDelete"
	}

	# Reset the collection information
	::cvis::queue::Add "$this ExecuteSetPreviousCollection"

	if {$DataIsVisible} {
	    ::cvis::queue::Add "$this ExecuteDataCreate"
	}
 	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}
    }
}


