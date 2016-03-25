##
## File:        cvisSlicePlane.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Slice plane
##

class cvisSlicePlane {
    inherit cvisDisplayObject

    # Slice plane
    variable SlicePlane "z"
    variable Slice 0

    # Data slice
    variable DataExists 0
    variable DataIsVisible 1

    variable TextureMap 1
    variable TextureMapInterpolate 0

    # Cell boundary
    variable CellBoundaryExists 0
    variable CellBoundaryIsVisible 0
    variable CellBoundaryWidth 1

    variable ScaleFactor 1.0
    variable TableRange "0.0 1.0"
    variable AxisScale "1.0 1.0 1.0"

    # Links to other objects
    variable Reader ""
    variable Interface ""
    variable LookupTable ""
    variable Depends ""
    variable PreviousCollection ""

    variable NumCellColors 8
    variable CellColors
    variable CellVisible

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	set NumCellColors 8
	# green
	set CellColors(0) "0 1 0"
	# yellow
	set CellColors(1) "1 1 0"
	# light blue
	set CellColors(2) "0 1 1"
	# purple
	set CellColors(3) "1 0 1"
	# pink 
	set CellColors(4) "1 0.75 0.80"
	# orange
	set CellColors(5) "1 0.65 0"
	# red
	set CellColors(6) "1 0 0"
	# blue
	set CellColors(7) "0 0 1"
	
	set CellVisible(0) "1"
	set CellVisible(1) "1"
	set CellVisible(2) "1"
	set CellVisible(3) "1"
	set CellVisible(4) "1"
	set CellVisible(5) "1"
	set CellVisible(6) "1"
	set CellVisible(7) "1"

    }

    method SetTextureMap {texturemap} {
	set tracervar [::itcl::local cvisTrace #auto]
	# SGS this should be schuduled and should be done
	# better?  Have one delete which tests for different
	# objects?
	if {$TextureMap != $texturemap} {
	    ExecuteDataDelete
	    set TextureMap $texturemap
	    ExecuteDataCreate
	}
    }

    method SetTextureMapInterpolate {texturemapinterpolate} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$TextureMapInterpolate != $texturemapinterpolate} {
	    set TextureMapInterpolate $texturemapinterpolate
	    if {$TextureMap && $DataExists} {
		set collection [$Reader GetCollection]
		$collection InitTraversal
		set item [$collection GetNextItem]
		set patch 0
		while {[string length $item]} { 
		    set item [$collection GetNextItem]
		}
	    }
	}
    }

    method SetPatchVisibility {level visible} {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method SetCellVisibility {level visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$CellVisible($level) != $visible} {
	    set CellVisible($level) $visible

	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    set patch 0
	    while {[string length $item]} {
		set patchlevel [$Reader GetPatchLevel $patch]
		
		if {$patchlevel == $level} {
		    if { [string length [info commands \
			    $this.$item.iedgeActor]] } {
			$this.$item.iedgeActor SetVisibility \
				[expr $CellVisible($level) \
				& $CellBoundaryIsVisible]
		    }
		}
		
		set item [$collection GetNextItem]
		incr patch
	    }
	}
    }

    method SetPatchColor {level r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method SetCellColor {level r g b} {
	set tracervar [::itcl::local cvisTrace #auto]

	set CellColors($level) "$r $g $b"

	if {[string length $Reader]} {

	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    set patch 0
	    while {[string length $item]} {
		set patchlevel [$Reader GetPatchLevel $patch]
		
		if {$patchlevel == $level} {
		    if { [string length [info commands \
			    $this.$item.iedgeActor]] } {
			eval [$this.$item.iedgeActor GetProperty] SetColor \
				$CellColors([expr $level % $NumCellColors])
		    }
		}
		
		set item [$collection GetNextItem]
		incr patch
	    }
	}
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
	ExecuteCellBoundarySetSlice
    }

    method SetSlicePlane {plane} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$SlicePlane == $plane} {
	    return
	}

	set SlicePlane $plane

	Update
    }

    method SetPatchBoundaryWidth {width} {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    #=======================================================================
    # CellBoundary 
    #=======================================================================
    method SetCellBoundaryWidth {width} {
	set tracervar [::itcl::local cvisTrace #auto]
	set CellBoundaryWidth $width

	# If no data reader then can't do anything
	if {![string length $Reader]} {
	    return
	}

	# If not visible then we don't need to do anything
	if {!$CellBoundaryIsVisible} {
	    return
	}

	if {!$CellBoundaryExists} {
	    return
	}

	set collection [[$Reader GetReader] GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    [$this.$item.iedgeActor GetProperty] SetLineWidth \
		    $CellBoundaryWidth
	    set item [$collection GetNextItem]
	}
    }

    method SetPatchBoundaryVisibility {vis} {
	set tracervar [::itcl::local cvisTrace #auto]
    }

    method SetCellBoundaryVisibility {vis} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	if {$vis == $CellBoundaryIsVisible} {
	    return
	}

	if {$CellBoundaryIsVisible} {
	    set CellBoundaryIsVisible 0
	    ::cvis::queue::Add "$this ExecuteCellBoundaryDelete"
	} {
	    set CellBoundaryIsVisible 1
	    ::cvis::queue::Add "$this ExecuteCellBoundaryCreate"
	}

    }

    method GetCellBoundaryVisibility {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $CellBoundaryIsVisible
    }

    method ExecuteCellBoundaryCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {![string length $RenderWindow] } { 
	    return
	}

	if {$CellBoundaryExists} {
	    return
	}

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {
	    scan [$item GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)

	    set level [$Reader GetPatchLevel $patch]
	    
	    
	    
	    

	    cvisLODActor $this.$item.iedgeActor
	    $this.$item.iedgeActor SetVisibility [expr $CellVisible($level) \
		    & $CellBoundaryIsVisible]
	    
	    set level [$Reader GetPatchLevel $patch]
	    
	    eval [$this.$item.iedgeActor GetProperty] SetColor \
		    $CellColors([expr $level % $NumCellColors])
	    
	    $RenderWindow AddActor $this.$item.iedgeActor

	    set item [$collection GetNextItem]
	    incr patch
	}

	set CellBoundaryExists 1

	$RenderWindow SetScaling

	ExecuteCellBoundarySetSlice
    }

    method ExecuteCellBoundaryDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set CellBoundaryExists 0	
	foreach item $PreviousCollection {
	    if {[string length \
		    [info commands $this.$item.iedgeActor]]} {
		
		$RenderWindow RemoveActor $this.$item.iedgeActor
		::itcl::delete object $this.$item.iedgeActor
	    }
	}
    }
    
    method ExecuteCellBoundarySetSlice {} {	
	set tracervar [::itcl::local cvisTrace #auto]

	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}

	# If no data reader then can't do anything
	if {![string length $Reader]} {
	    return
	}

	# If not visible then we don't need to do anything
	if {!$CellBoundaryIsVisible} {
	    return
	}

	if {!$CellBoundaryExists} {
	    return
	}

	scan [[$Reader GetReader] GetBounds] "%f %f %f %f %f %f" \
		real(min,x) real(max,x) \
		real(min,y) real(max,y) \
		real(min,z) real(max,z)

	# Push the slice in from the boundaries
	# SGS is this really needed anymore for cell centered data?
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
	set patch 0
	while {[string length $item]} {
	    set level [$Reader GetPatchLevel $patch]

	    scan [$item GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)
	    
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    item_bounds(min,x) item_bounds(max,x) \
		    item_bounds(min,y) item_bounds(max,y) \
		    item_bounds(min,z) item_bounds(max,z)

	    if { $Slice >= $item_bounds(min,$SlicePlane) && \
		    $Slice <= $item_bounds(max,$SlicePlane) } {

		switch $SlicePlane {
		    "x" {
				1 $dim(y) $dim(z) 
				$Slice $item_bounds(min,y) $item_bounds(min,z)
		    } 
		    "y" {
				$dim(x) 1 $dim(z) 
				$item_bounds(min,x) $Slice $item_bounds(min,z)
		    } 
		    "z" {
				$dim(x) $dim(y) 1 
				$item_bounds(min,x) $item_bounds(min,y) $Slice
		    } 
		}
		
		$this.$item.iedgeActor SetVisibility  \
			[expr $CellVisible($level) & $CellBoundaryIsVisible]
	    } {
		$this.$item.iedgeActor SetVisibility 0
	    }

	    [$this.$item.iedgeActor GetProperty] SetLineWidth \
		    $CellBoundaryWidth

	    set item [$collection GetNextItem]
	    incr patch
	}
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
	if {$TextureMap} { 
	    ExecuteTextureDataCreate
	} {
	    ExecuteTriangleDataCreate
	}
    }

    method ExecuteDataDelete {} {
	if {$TextureMap} { 
	    ExecuteTextureDataDelete
	} {
	    ExecuteTriangleDataDelete
	}
    }

    method ExecuteDataSetSlice {} {
	if {$TextureMap} { 
	    ExecuteTextureDataSetSlice
	} {
	    ExecuteTriangleDataSetSlice
	}
    }

    #=======================================================================
    # Triangle Data Slice
    #=======================================================================
    method ExecuteTriangleDataCreate {} {
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
	    
		    GetCenter]
	    
		    $item
	    
	    

		    $::cvis::options::LOWMEM
	    
	    if {[string length $LookupTable] } { 
			[$LookupTable GetLookupTable]
	    }
	    
	    cvisLODActor $this.$item.planeActor
	    
	    # Turn Actor off by default since it might not lie on the
	    # current slice plane.
	    $this.$item.planeActor SetVisibility 0

	    $RenderWindow AddActor $this.$item.planeActor
	    
	    set item [$collection GetNextItem]
	}
	
	$RenderWindow SetScaling
	
	ExecuteDataSetSlice 
    }

    method ExecuteTriangleDataDelete { } {
	set tracervar [::itcl::local cvisTrace #auto]

	set DataExists 0
	foreach item $PreviousCollection {
	    if {[string length \
		
		$RenderWindow RemoveActor $this.$item.planeActor
		::itcl::delete object $this.$item.planeActor
	    }
	}
    }

    method ExecuteTriangleDataSetSlice {} {
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
	# SGS is this really needed anymore for cell centered data?
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
		
		
		switch $SlicePlane {
		    "x" {
		    }
		    "y" {
		    }
		    "z" {
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

    method ExecuteTextureDataCreate {} {
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



	    # There is a problem with VTK texture mapping if one of
	    # the dims is 1 since VTK uses dim=1 to indicate the 
	    # dimension which the 2D image is lying in.  So get 
	    # around this by replicating the data.
	    # 
	    if { ([string compare $SlicePlane "z"] == 0  || \
		    [string compare $SlicePlane "y"] == 0) && $dim(x) == 2} {
	    } {
	    }


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

	    # texture from reslice filter
	    


#	    global LOWMEM
	    
	    if {[string length $LookupTable] } { 
			[$LookupTable GetLookupTable]
	    }
	    
	    cvisActor $this.$item.planeActor

	    # Turn Actor off by default since it might not lie on the
	    # current slice plane.
	    $this.$item.planeActor SetVisibility 0
	    
	    $RenderWindow AddActor $this.$item.planeActor

	    set item [$collection GetNextItem]
	}
	
	$RenderWindow SetScaling
	
	ExecuteDataSetSlice 
    }


    method ExecuteTextureDataDelete { } {
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

    method ExecuteTextureDataSetSlice {} {
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
	# SGS is this really needed anymore for cell centered data?
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
    
    method SetScaleFactor {factor} {
	set tracervar [::itcl::local cvisTrace #auto]
	set ScaleFactor $factor

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    set item [$collection GetNextItem]
	}
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
	    set Slice [$Interface GetSliceInScreenCoord]
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	# Schedule CPU intensive tasks that accomplish the update
	if {$CellBoundaryExists} {
	    set CellBoundaryExists 0
	    ::cvis::queue::Add "$this ExecuteCellBoundaryDelete"
	}
	if {$DataExists} {
	    set DataExists 0
	    ::cvis::queue::Add "$this ExecuteDataDelete"
	}

	# Reset the collection information
	::cvis::queue::Add "$this ExecuteSetPreviousCollection"

	if {$DataIsVisible} {
	    ::cvis::queue::Add "$this ExecuteDataCreate"
	}
	if {$CellBoundaryIsVisible} {
	    ::cvis::queue::Add "$this ExecuteCellBoundaryCreate"
	}

 	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}
    }
}


