#
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisVector.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Slice plane
##

class cvisVector {
    inherit cvisDisplayObject

    # Slice plane
    variable VOI

    # Data slice
    variable DataExists 0
    variable DataIsVisible 1

    variable ScalarVisibility 1
    variable VectorColor "1 1 1"

    variable ScaleFactor 1.0
    variable TableRange "0.0 1.0"
    variable AxisScale "1.0 1.0 1.0"

    # Links to other objects
    variable Reader ""
    variable Interface ""
    variable LookupTable ""
    variable Depends ""
    variable PreviousCollection ""

    variable ScalarRange

    variable Outputs ""

    variable Crop 0
    variable CroppingRange 
    variable VectorScaleFactor

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

    method SetVectorColor { r g b} {
	set tracervar [::itcl::local cvisTrace #auto]

	set VectorColor "$r $g $b"
	
	if {[string length $Reader]} {
	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		if { [string length [info commands \
			$this.$item.planeActor]] } {
		    eval [$this.$item.planeActor GetProperty] SetColor $VectorColor
		}
		set item [$collection GetNextItem]
	    }
	}
    }

    method GetVectorColor {} {
	return $VectorColor
    }

    method GetScalarVisibilty {} {
	return $ScalarVisibility
    }

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	set ScalarRange(min) 0
	set ScalarRange(max) 1

	set CroppingRange "0 1"
	set VectorScaleFactor "0 1"
    }

    #=======================================================================
    # Access Functions
    #=======================================================================
    method GetScalarRange {} {
	return "$ScalarRange(min) $ScalarRange(max)"
    }

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

    method SetVOI {min_x max_x min_y max_y min_z max_z} {
	set VOI(x,min) $min_x
	set VOI(x,max) $max_x
	set VOI(y,min) $min_y
	set VOI(y,max) $max_y
	set VOI(z,min) $min_z
	set VOI(z,max) $max_z

	::cvis::queue::Add "$this ExecuteSetVOI"
    }

    method ExecuteSetVOI {} {
	ExecuteDataSetSlice
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
	    ExecuteDataDelete
	} {
	    set DataIsVisible 1
	    ExecuteDataCreate
	}
    }
    
    method GetDataVisibility {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $DataIsVisible
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

#***********************************************************************

    method ExecuteDataCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	if {$DataExists} {
	    return
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	set ScalarRange(min) $::cvis::math::FLT_MAX
	set ScalarRange(max) -$::cvis::math::FLT_MAX

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

		    $origin_new(x) \
		    $origin_new(y) \
		    $origin_new(z)
		    
		    $spacing(x) \
		    $spacing(y) \
		    $spacing(z) 



		    min max

	    if {$min < $ScalarRange(min)} {
		set ScalarRange(min) $min
	    }

	    if {$max > $ScalarRange(max)} {
		set ScalarRange(max) $max
	    }






	    if {[string length $LookupTable] } { 
			[$LookupTable GetLookupTable]
	    }

	    cvisActor $this.$item.planeActor

	    # Turn Actor off by default since it might not lie on the
	    # current slice plane.
	    $this.$item.planeActor SetVisibility 0

	    eval [$this.$item.planeActor GetProperty] SetColor $VectorColor
	    
	    $RenderWindow AddActor $this.$item.planeActor

	    set item [$collection GetNextItem]
	}

	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    set item [$collection GetNextItem]
	}

	$RenderWindow SetScaling
	
	ExecuteDataSetSlice 

	# Update children of this module since Scalar/MinMax may have
	# Changed.
	foreach i $Outputs { 
	    $i SetInitialize
	    $i Update
	}
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

	    foreach plane {x y z} i {0 1 2}  {
		set origin($plane) [lindex [$item GetOrigin] $i]
		set dim($plane)    [lindex [$item GetDimensions] $i]
		set spacing($plane)    [lindex [$item GetSpacing] $i]
	    }

	    if { ($VOI(x,min) > $real(max,x) || \
		    $VOI(x,max) < $real(min,x)) || \
		    ($VOI(y,min) > $real(max,y) || \
		    $VOI(y,max) < $real(min,y)) || \
		    ($VOI(z,min) > $real(max,z) || \
		    $VOI(z,max) < $real(min,z)) } {
		$this.$item.planeActor SetVisibility 0
	    } {
		foreach plane {x y z} {
		    set IndexVOI($plane,min) [expr int(($VOI($plane,min)-$origin($plane)) / \
			    $spacing($plane))]

		    if { $IndexVOI($plane,min) < 0 } { 
			set IndexVOI($plane,min) 0
		    }
			
		    set IndexVOI($plane,max) [expr int(($VOI($plane,max)-$origin($plane)) / \
			    $spacing($plane))]

		    if { $IndexVOI($plane,max) > $dim($plane) } { 
			set IndexVOI($plane,max) $dim($plane)
		    }
		}

				$IndexVOI(x,min) $IndexVOI(x,max) \
				$IndexVOI(y,min) $IndexVOI(y,max) \
				$IndexVOI(z,min) $IndexVOI(z,max) 

		$this.$item.planeActor SetVisibility $DataIsVisible
		$this.$item.planeActor Update
	    }
	    set item [$collection GetNextItem]
	}
    }

#***********************************************************************

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
    
    # SGS new stuff
    method GetCrop {} {
	return $Crop
    }

    method SetCrop {flag} {
	set Crop $flag
	::cvis::queue::Add "$this ExecuteSetCrop"
    }

    method ExecuteSetCrop {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	    while {[string length $item]} {
		set item [$collection GetNextItem]
	    }
	}
    }

    method GetCroppingRange {} {
	return $CroppingRange
    }

    method SetCroppingRange {min max} {
	set CroppingRange "$min $max"
	::cvis::queue::Add "$this ExecuteSetCroppingRange"
    }

    method ExecuteSetCroppingRange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	    while {[string length $item]} {
		set item [$collection GetNextItem]
	    }

	}
    }

    method GetVectorScaleFactor {} {
	return $VectorScaleFactor
    }

    method SetVectorScaleFactor {min max} {
	set VectorScaleFactor "$min $max"
	::cvis::queue::Add "$this ExecuteSetVectorScaleFactor"
    }

    method ExecuteSetVectorScaleFactor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	    while {[string length $item]} {
		set item [$collection GetNextItem]
	    }
	}
    }

    # End SGS new stuff

    method GetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $TableRange 
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

	if {$DataExists} {
	    set DataExists 0
	    ExecuteDataDelete
	}
    
	# Reset the collection information
	ExecuteSetPreviousCollection
    
	if {$DataIsVisible} {
	    ExecuteDataCreate
	}

	# Update user interface elements
	if {[string length $Interface]} {
	    $Interface Update
	    # Get values that might have changed based on UI control
	}

	if {![string length $RenderWindow] } { 
	    return
	}

#
# SGS this needs to be reworked, we need the scalar range to
# update the user interface so we must do this before 
# the update.....on the other hand we would like to schedule
# this...what to do?
#
#  	if {0} {
#  	    # Schedule CPU intensive tasks that accomplish the update
#  	    if {$CellBoundaryExists} {
#  		set CellBoundaryExists 0
#  		::cvis::queue::Add "$this ExecuteCellBoundaryDelete"
#  	    }
#  	    if {$DataExists} {
#  		set DataExists 0
#  		::cvis::queue::Add "$this ExecuteDataDelete"
#  	    }
	    
#  	    # Reset the collection information
#  	    ::cvis::queue::Add "$this ExecuteSetPreviousCollection"
	    
#  	    if {$DataIsVisible} {
#  		::cvis::queue::Add "$this ExecuteDataCreate"
#  	    }
#  	    if {$CellBoundaryIsVisible} {
#  		::cvis::queue::Add "$this ExecuteCellBoundaryCreate"
#  	    }
#  	}
	
 	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}

	# Update children of this module
	foreach i $Outputs { 
	    $i SetInitialize
	    $i Update
	}
    }

    method AddOutput {list} {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Outputs $list
    }
}

