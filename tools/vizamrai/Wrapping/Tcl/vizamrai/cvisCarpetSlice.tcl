##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisCarpetSlice.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Carpet Slice
##

class cvisCarpetSlice {
    inherit cvisDisplayObject

    variable DataIsVisible 1
    
    variable Reader ""
    variable SlicePlane "z"
    variable Slice 0.5

    variable PreviousCollection ""
    variable PreviousLevel 0

    variable Offset 0.0
    variable ScaleFactor 0.5
    variable NormalizedScaleFactor 1.0
    variable TableRange "0.0 1.0"

    variable LookupTable ""

    variable Interface ""

    variable AxisScale "1.0 1.0 1.0"

    variable BoundaryEdgesVisible 1
    variable BoundaryEdgeWidth 1
    variable BoundaryEdgesInterior 1
    variable BoundaryExists 0

    variable LevelHasData

    variable CellToNode 0
    
    variable DataBrowser ""

    variable BoundingBoxes ""

    variable AddPatches 0

    variable DataExists 0

    variable NumBorderColors 8
    variable BorderColors
    variable BorderVisible

    # Flag to help with flow.  On a new file load we are getting two
    # inputs, lookup table and the reader.  Don't want to 
    # set scale values for the old data since we are going to delete it
    # next.
    variable DataCreateMode 0

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	if {$AddPatches} {
	}

	set NumBorderColors 8
	# green
	set BorderColors(0) "0 1 0"
	# yellow
	set BorderColors(1) "1 1 0"
	# light blue
	set BorderColors(2) "0 1 1"
	# purple
	set BorderColors(3) "1 0 1"
	# pink 
	set BorderColors(4) "1 0.75 0.80"
	# orange
	set BorderColors(5) "1 0.65 0"
	# red
	set BorderColors(6) "1 0 0"
	# blue
	set BorderColors(7) "0 0 1"
	
	set BorderVisible(0) "1"
	set BorderVisible(1) "1"
	set BorderVisible(2) "1"
	set BorderVisible(3) "1"
	set BorderVisible(4) "1"
	set BorderVisible(5) "1"
	set BorderVisible(6) "1"
	set BorderVisible(7) "1"
    }

    method SetBorderVisibility {level visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BorderVisible($level) != $visible} {
	    set BorderVisible($level) $visible
	    if {[string length $Reader]} {
		if {$BoundaryEdgesInterior} {
		    set collection [$Reader GetCollection]
		    $collection InitTraversal
		    set patch 0
		    set item [$collection GetNextItem]
		    while {[string length $item]} {
			if {[$Reader GetPatchLevel $patch] == $level} {
			    if { [string length [info commands \
				    $this.$item.boundaryActor]] } {
				eval [$this.$item.boundaryActor SetVisibility \
					[expr $BoundaryEdgesVisible & \
					$BorderVisible([expr $level % \
					$NumBorderColors])]]
			    }
			}
			incr patch
			set item [$collection GetNextItem]
		    }
	        } {
		    if {$LevelHasData($level)} {
			if { [string length [info commands \
				$this.$level.boundaryActor]]} {
			    $this.$level.boundaryActor SetVisibility \
				    [expr $BoundaryEdgesVisible & \
				    $BorderVisible([expr $level % \
					$NumBorderColors])]
			}
		    } 
		}
	    }
	}
    }

    method SetBorderColor {level r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BorderColors($level) != "$r $g $b"} {
	    set BorderColors($level) "$r $g $b"
	    if {[string length $Reader]} {
		if {$BoundaryEdgesInterior} {
		    set collection [$Reader GetCollection]
		    $collection InitTraversal
		    set patch 0
		    set item [$collection GetNextItem]
		    while {[string length $item]} {
			if {[$Reader GetPatchLevel $patch] == $level} {
			    if { [string length [info commands \
				    $this.$item.boundaryActor]] } {
				
				eval [$this.$item.boundaryActor GetProperty] \
					SetColor $BorderColors([expr $level % \
					$NumBorderColors])
			    }
			}
			incr patch
			set item [$collection GetNextItem]
		    }
		} {
		    if {$LevelHasData($level)} {
			if {[string length [info commands \
				$this.$level.boundaryActor]]} {
			    eval [$this.$level.boundaryActor GetProperty] \
				    SetColor $BorderColors([expr $level % \
				    $NumBorderColors])
			}
		    }
		}
	    }
	}
    }

    method SetDataBrowser {browser} {
	set tracervar [::itcl::local cvisTrace #auto]
	set DataBrowser $browser
    }

    method GetDataBrowser { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $DataBrowser
    }

    method GetScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [$Reader GetScaling]
	} {
	    return 1.0
	}
    }

    method SetCellToNode {value} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$CellToNode != $value} {
	    set CellToNode $value
	    Update
	}
    }
    
    method GetCellToNode {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $CellToNode
    }
    
    method WriteCut {filename} {
	set tracervar [::itcl::local cvisTrace #auto]

	vtkPolyDataWriter writer
	writer SetFileTypeToASCII

	set collection [[$Reader GetReader] GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {
	    writer SetFileName $filename.$patch
	    writer Write
	    incr patch
	    set item [$collection GetNextItem]
	}

	writer Delete
    }

    method SetBoundaryEdgesVisibility {visibility} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BoundaryEdgesVisible != $visibility} {
	    set BoundaryEdgesVisible $visibility
	    if {[string length $Reader]} {
		if {$BoundaryEdgesInterior} {
		    set collection [$Reader GetCollection]
		    $collection InitTraversal
		    set patch 0
		    set item [$collection GetNextItem]
		    while {[string length $item]} {
			set level [$Reader GetPatchLevel $patch]
			if { [string length [info commands \
				$this.$item.boundaryActor]] } {
			    eval [$this.$item.boundaryActor SetVisibility \
				    [expr $BoundaryEdgesVisible \
				    & $BorderVisible([expr $level % \
					$NumBorderColors])]]
			}
			incr patch
			set item [$collection GetNextItem]
		    }
		    } {
		    for {set level 0} \
			    {$level < [[$Reader GetReader] GetNumberOfLevels]} \
			    {incr level} {
			if {$LevelHasData($level)} {
			    if {[string length [info commands \
				    $this.$level.boundaryActor]]} {
				$this.$level.boundaryActor SetVisibility \
					[expr $BoundaryEdgesVisible & \
					$BorderVisible([expr $level % \
					$NumBorderColors])]
			    }
			}
		    }
		}
	    }
	}
    }

    method GetBoundaryEdgesVisible { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $BoundaryEdgesVisible
    }

    method SetBoundaryEdgesInterior {value} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BoundaryEdgesInterior == $value} {
	    return
	}

	if {$BoundaryEdgesInterior} {
	    set BoundaryEdgesInterior 0
	    ::cvis::queue::Add "$this ExecuteBoundaryDelete"
	    ::cvis::queue::Add "$this ExecuteBoundaryCreate"
	} {
	    set BoundaryEdgesInterior 1
	    ::cvis::queue::Add "$this ExecuteBoundaryDelete"
	    ::cvis::queue::Add "$this ExecuteBoundaryCreate"
	}
    }

    method GetBoundaryEdgesInterior { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $BoundaryEdgesInterior
    }

    
    method SetBoundaryEdgeWidth { width } {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$BoundaryEdgeWidth == $width} {
	    return
	}

	set BoundaryEdgeWidth $width

	if {[string length $Reader]} {
	    ::cvis::queue::Add "$this ExecuteSetBoundaryEdgeWidth"
	}
    }

    method ExecuteSetBoundaryEdgeWidth {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$BoundaryExists} {
	    if {$BoundaryEdgesInterior} {
		set collection [[$Reader GetReader] GetCollection]
		$collection InitTraversal
		set item [$collection GetNextItem]
		while {[string length $item]} {
		    if { [string length \
			[$this.$item.boundaryActor GetProperty] \
				SetLineWidth $BoundaryEdgeWidth 
		    }
		    set item [$collection GetNextItem]
		}
	    } {
		for {set level 0} \
			{$level < [[$Reader GetReader] GetNumberOfLevels]} \
			{incr level} {
		    if {$LevelHasData($level)} {
			[$this.$level.boundaryActor GetProperty] \
				SetLineWidth $BoundaryEdgeWidth 
		    }
		}
	    }
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
	if {$AddPatches} {
	}
    }

    method SetInterface { interface } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Interface $interface
	$interface Update
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

    method GetSlice {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Slice
    }

    method ExecuteDataDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set DataExists 0

	foreach item $PreviousCollection {
	    if {[string length \
		
		$RenderWindow RemoveActor $this.$item.planeActor
		::itcl::delete object $this.$item.planeActor 
	    }

	    if {[string length \
	    }

	    if {[string length \
	    }

	    if {[string length \
	    }
	}
    }

    method ExecuteBoundaryDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set BoundaryExists 0
	for {set level 0} {$level < $PreviousLevel} {incr level} {
		$RenderWindow RemoveActor $this.$level.boundaryActor
		::itcl::delete object $this.$level.boundaryActor 
	    }
	    
	    }
	}
	
	foreach item $PreviousCollection {
	    if {[string length \
		    [info commands $this.$item.boundaryActor]]} {
		$RenderWindow RemoveActor $this.$item.boundaryActor
		::itcl::delete object $this.$item.boundaryActor
	    }
	    
	    if {[string length \
	    }
	    
	    if {[string length \
		
	    }
	}
    }
    
    method ExecuteBoundaryCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$BoundaryEdgesInterior} {
	    ExecuteInteriorBoundaryCreate
	} {
	    ExecuteAllBoundaryCreate
	}
    }

    method ExecuteAllBoundaryCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {![string length $RenderWindow] } { 
	    return
	}

	if {$BoundaryExists} {
	    return
	}

	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}

	set BoundaryExists 1
	for {set level 0} {$level < [[$Reader GetReader] GetNumberOfLevels]} \
		{incr level} {
	    set LevelHasData($level) 0
	}

	scan [[$Reader GetReader] GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z) 

	set offset [expr $Offset - [lindex $TableRange 0]]
	switch $SlicePlane {
	    "x" {
		set scale_factor   "$ScaleFactor 1 1"
		set actor_position "[expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0 0"
		set actor_origin   "$Slice 0 0"
	    }
	    "y" {
		set scale_factor   "1 $ScaleFactor 1"
		set actor_position "0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0"
		set actor_origin   "0 $Slice 0"
	    }
	    "z" {
		set scale_factor   "1 1 $ScaleFactor"
		set actor_position "0 0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)]"
		set actor_origin   "0 0 $Slice"
	    }
	}

	set Scaling [[$Reader GetReader] GetScaling] 

	# Set to value which is offset if looking at the ends of
	# the box othewise we might not intersect
	set jitter [expr $IndexScaling($SlicePlane) * $Scaling * \
		$::cvis::math::BOUNDARY_JITTER ]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {
	    set level [$Reader GetPatchLevel $patch]
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    real(min,x) real(max,x) \
		    real(min,y) real(max,y) \
		    real(min,z) real(max,z)

	    scan [$item GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)
	    
	    scan [$item GetSpacing] "%f %f %f" \
		    spacing(x) spacing(y)  spacing(z)
	    
	    if {$CellToNode} {
		set intersect_boundary(max) $real(max,$SlicePlane)
		set intersect_boundary(min) $real(min,$SlicePlane)
	    } {
		if {$dim($SlicePlane) == 2} {
		    set intersect_boundary(max) $real(max,$SlicePlane)
		    set intersect_boundary(min) $real(min,$SlicePlane)
		} {
		    set intersect_boundary(max) [expr $real(max,$SlicePlane) \
			    - $spacing($SlicePlane) /2.0]
		    set intersect_boundary(min) [expr $real(min,$SlicePlane) \
			    + $spacing($SlicePlane) / 2.0]
		}
	    }

	    if { $Slice >= [expr $intersect_boundary(min) - $jitter] && \
		    $Slice <= [expr $intersect_boundary(max) + $jitter]} {

		if { $Slice < [expr $intersect_boundary(min) + $jitter] } {
		    set Slice [expr $Slice + $jitter]
		}
		
		if { $Slice > [expr $intersect_boundary(max) - $jitter] } {
		    set Slice [expr $Slice - $jitter]
		}
		
		set level [$Reader GetPatchLevel $patch]

		# Construct Exterior boundary edges for this
		# patch
		
		
		
		
		foreach plane "x y z" {
		    set rotate($plane) 0.0
		}
		
		if { $SlicePlane == "x" } {
		    set rotate(y) 90.0
		}
		
		if { $SlicePlane == "y" } {
		    set rotate(x) 90.0
		}
		
		
		
		
		
		if { $SlicePlane == "x" } {
		    set rotate(y) -90.0
		}
		
		if { $SlicePlane == "y" } {
		    set rotate(x) -90.0
		}
		
		
		
		
		
		set LevelHasData($level) 1
		
	    }

	    set item [$collection GetNextItem]
	    incr patch
	}

	set PreviousLevel [[$Reader GetReader] GetNumberOfLevels]
	for {set level 0} \
		{$level < [[$Reader GetReader] GetNumberOfLevels]} \
		{incr level} {
	    if {$LevelHasData($level)} {
			$::cvis::math::BOUNDARY_TOLERANCE
		
		
		
			$::cvis::options::LOWMEM

		
		cvisWarpActor $this.$level.boundaryActor
		eval [$this.$level.boundaryActor GetProperty] SetColor \
			$BorderColors([expr $level % $NumBorderColors])
		[$this.$level.boundaryActor GetProperty] SetAmbient  1
		[$this.$level.boundaryActor GetProperty] SetDiffuse  0
		[$this.$level.boundaryActor GetProperty] SetSpecular 0
		[$this.$level.boundaryActor GetProperty] SetLineWidth \
			$BoundaryEdgeWidth

		eval $this.$level.boundaryActor SetScaleFactor $scale_factor

		$this.$level.boundaryActor SetVisibility \
			[expr $BoundaryEdgesVisible & \
			$BorderVisible([expr $level % \
					$NumBorderColors])]

		eval $this.$level.boundaryActor SetPosition $actor_position
		eval $this.$level.boundaryActor SetOrigin   $actor_origin
		
		$this.$level.boundaryActor SetMapper \
		
		$RenderWindow AddActor $this.$level.boundaryActor
	    }
	}
    }

    method ExecuteInteriorBoundaryCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {![string length $RenderWindow] } { 
	    return
	}

	if {$BoundaryExists} {
	    return
	}

	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}


	set BoundaryExists 1

	scan [[$Reader GetReader] GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z) 

	set offset [expr $Offset - [lindex $TableRange 0]]
	switch $SlicePlane {
	    "x" {
		set scale_factor "$ScaleFactor 1 1"
		set actor_position "[expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0 0"
		set actor_origin   "$Slice 0 0"
	    }
	    "y" {
		set scale_factor "1 $ScaleFactor 1"
		set actor_position "0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0"
		set actor_origin   "0 $Slice 0"
	    }
	    "z" {
		set scale_factor "1 1 $ScaleFactor"
		set actor_position "0 0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)]"
		set actor_origin   "0 0 $Slice"
	    }
	}


	set Scaling [[$Reader GetReader] GetScaling] 

	# Set to value which is offset if looking at the ends of
	# the box othewise we might not intersect
	set jitter [expr $IndexScaling($SlicePlane) * $Scaling * \
		$::cvis::math::BOUNDARY_JITTER ]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {
	    set level [$Reader GetPatchLevel $patch]
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    real(min,x) real(max,x) \
		    real(min,y) real(max,y) \
		    real(min,z) real(max,z)

	    scan [$item GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)
	    
	    scan [$item GetSpacing] "%f %f %f" \
		    spacing(x) spacing(y)  spacing(z)
	    
	    if {$CellToNode} {
		set intersect_boundary(max) $real(max,$SlicePlane)
		set intersect_boundary(min) $real(min,$SlicePlane)
	    } {
		if {$dim($SlicePlane) == 2} {
		    set intersect_boundary(max) $real(max,$SlicePlane)
		    set intersect_boundary(min) $real(min,$SlicePlane)
		} {
		    set intersect_boundary(max) [expr $real(max,$SlicePlane) \
			    - $spacing($SlicePlane) /2.0]
		    set intersect_boundary(min) [expr $real(min,$SlicePlane) \
			    + $spacing($SlicePlane) / 2.0]
		}
	    }

	    if { $Slice >= [expr $intersect_boundary(min) - $jitter] && \
		    $Slice <= [expr $intersect_boundary(max) + $jitter]} {

		if { $Slice < [expr $intersect_boundary(min) + $jitter] } {
		    set Slice [expr $Slice + $jitter]
		}
		
		if { $Slice > [expr $intersect_boundary(max) - $jitter] } {
		    set Slice [expr $Slice - $jitter]
		}

		# Boundary Outline
		
			$::cvis::options::LOWMEM
		cvisWarpActor $this.$item.boundaryActor
		
		eval [$this.$item.boundaryActor GetProperty] SetColor \
			$BorderColors([expr $level % $NumBorderColors])
		[$this.$item.boundaryActor GetProperty] SetAmbient  1
		[$this.$item.boundaryActor GetProperty] SetDiffuse  0
		[$this.$item.boundaryActor GetProperty] SetSpecular 0
		[$this.$item.boundaryActor GetProperty] SetLineWidth \
			$BoundaryEdgeWidth

		eval $this.$item.boundaryActor SetScaleFactor $scale_factor
		
		$this.$item.boundaryActor SetVisibility \
			[expr $BoundaryEdgesVisible && \
			$BorderVisible([expr $level % \
					$NumBorderColors])]

		eval $this.$item.boundaryActor SetPosition $actor_position
		eval $this.$item.boundaryActor SetOrigin   $actor_origin

		$this.$item.boundaryActor SetMapper \

		$RenderWindow AddActor $this.$item.boundaryActor
	    } 
	    
	    set item [$collection GetNextItem]
	    incr patch
	}
    }

    method ExecuteDataCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {![string length $RenderWindow] } { 
	    return
	}

	if {$DataExists} {
	    return
	}

	if {$AddPatches} {
	}

	set DataExists 1
	
	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}

	ComputeNormalizedScaleFactor
	
	scan [[$Reader GetReader] GetBounds] "%f %f %f %f %f %f" \
		real(min,x) real(max,x) \
		real(min,y) real(max,y) \
		real(min,z) real(max,z)

	scan [[$Reader GetReader] GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z) 

	set offset [expr $Offset - [lindex $TableRange 0]]
	switch $SlicePlane {
	    "x" {
		set scale_factor "$ScaleFactor 1 1"
		set actor_position "[expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0 0"
		set actor_origin   "$Slice 0 0"
	    }
	    "y" {
		set scale_factor "1 $ScaleFactor 1"
		set actor_position "0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0"
		set actor_origin   "0 $Slice 0"
	    }
	    "z" {
		set scale_factor "1 1 $ScaleFactor"
		set actor_position "0 0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)]"
		set actor_origin   "0 0 $Slice"
	    }
	}

	set Scaling [[$Reader GetReader] GetScaling] 

	# Set to value which is offset if looking at the ends of
	# the box othewise we might not intersect
	set jitter [expr $IndexScaling($SlicePlane) * $Scaling * \
		$::cvis::math::BOUNDARY_JITTER ]

	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	set patch 0
	while {[string length $item]} {
	    set level [$Reader GetPatchLevel $patch]
	    
	    lappend PreviousCollection $item

	    # Need to do more complicated stuff here for 
	    # see the vizamrai for more information
	    # need to set to ends and set vizibility off
	    
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    real(min,x) real(max,x) \
		    real(min,y) real(max,y) \
		    real(min,z) real(max,z)

	    scan [$item GetDimensions] "%d %d %d" \
		    dim(x) dim(y) dim(z)
	    
	    scan [$item GetSpacing] "%f %f %f" \
		    spacing(x) spacing(y)  spacing(z)
	    
	    scan [$item GetOrigin] "%f %f %f" \
		    origin(x) origin(y)  origin(z)
	    
	    # SGS Changed this but is this a bug or not?
	    # should work on intersection routines
	    # was test on if {$CellToNode} which caused datasets
	    # with dim=1 to fail since min==max when moved half cell
	    # in
	    if {$CellToNode} {
		set intersect_boundary(max) $real(max,$SlicePlane)
		set intersect_boundary(min) $real(min,$SlicePlane)
	    } {
		if {$dim($SlicePlane) == 2} {
		    set intersect_boundary(max) $real(max,$SlicePlane)
		    set intersect_boundary(min) $real(min,$SlicePlane)
		} {
		    set intersect_boundary(max) [expr $real(max,$SlicePlane) \
			    - $spacing($SlicePlane) /2.0]
		    set intersect_boundary(min) [expr $real(min,$SlicePlane) \
			    + $spacing($SlicePlane) / 2.0]
		}
	    }

	    if { $Slice >= [expr $intersect_boundary(min) - $jitter] && \
		    $Slice <= [expr $intersect_boundary(max) + $jitter]} {

		if { $Slice < [expr $intersect_boundary(min) + $jitter] } {
		    set Slice [expr $Slice + $jitter]
		}
		
		if { $Slice > [expr $intersect_boundary(max) - $jitter] } {
		    set Slice [expr $Slice - $jitter]
		}

		if {$CellToNode} {
		    
		    switch $SlicePlane {
			"x" {
			}
			"y" {
			}
			"z" {
			}
		    }
		    
		    

		} {

		    $item Update

		    # Convert to Cell to Node centered data

		    # Shrink size
		    # if one D can bypass the cut and display
		    # if it is in the same plane otherwise we get a line
		    foreach plane "x y z" {
			if { $dim($plane) > 1 } {
			    incr dim($plane) -1
			}

			set cell_centered_origin($plane) \
				[expr $origin($plane) + \
				$spacing($plane)/2.0]
		    }

			    $spacing(x) $spacing(y) $spacing(z)

			    $cell_centered_origin(x) \
			    $cell_centered_origin(y) \
			    $cell_centered_origin(z)
		    
			    [[$item GetCellData] GetScalars]


		    # Grab the plane
		    if { $dim($SlicePlane) == 1 } {

		    } {


			switch $SlicePlane {
			    "x" {
			    } 
			    "y" {
			    }
			    "z" {
			    }
			}



			
		    }
		}

		# Warp it
		
		
			$::cvis::options::LOWMEM
		
		
		if {[string length $LookupTable] } { 
			    [$LookupTable GetLookupTable]
		}

		cvisLODWarpActor $this.$item.planeActor
		
		eval $this.$item.planeActor SetScaleFactor $scale_factor

		eval $this.$item.planeActor SetPosition $actor_position
		eval $this.$item.planeActor SetOrigin   $actor_origin

		$RenderWindow AddActor $this.$item.planeActor
	    }
	    
	    set item [$collection GetNextItem]
	    incr patch
	}


	#-----------------------------------------------------
	$RenderWindow SetScaling

	set DataCreateMode 0
    }

     method SetSlice {value} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	if {[expr $value != $Slice]} {
	    set Slice $value

	    Update

	    if {$AddPatches} {
	    }
	}
    }

    method SetSlicePlane {plane} {
	set tracervar [::itcl::local cvisTrace #auto]
	set SlicePlane $plane

	Update
    }

    method SetRenderWindow {window} {
	set tracervar [::itcl::local cvisTrace #auto]
	set RenderWindow $window
	if {$AddPatches} {
	}
    }
    
    method SetOffset {offset} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$Offset != $offset} {
	    set Offset $offset
	    
	    if [string length $Reader] { 

		set offset [expr $Offset - [lindex $TableRange 0]]
		switch $SlicePlane {
		    "x" {
			set actor_position "[expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0 0"
		    }
		    "y" {
			set actor_position "0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0"
		    }
		    "z" {
			set actor_position "0 0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)]"
		    }
		}
	
		set collection [$Reader GetCollection]
		$collection InitTraversal
		set item [$collection GetNextItem]
		while {[string length $item]} {
		    if { [string length [info commands $this.$item.planeActor]] } {
			eval $this.$item.planeActor SetPosition $actor_position
		    }
		    set item [$collection GetNextItem]
		}
		
		
		if {$BoundaryEdgesInterior} {
		    set collection [$Reader GetCollection]
		    $collection InitTraversal
		    set patch 0
		    set item [$collection GetNextItem]
		    while {[string length $item]} {
			if { [string length [info commands \
			    eval $this.$item.boundaryActor SetPosition $actor_position
			}
			incr patch
			set item [$collection GetNextItem]
		    }
		} {
		    for {set level 0} \
			    {$level < [[$Reader GetReader] GetNumberOfLevels]} \
			    {incr level} {
			if {$LevelHasData($level)} {
			    if {[string length [info commands \
				eval $this.$level.boundaryActor SetPosition $actor_position
			    }
			}
		    }
		}
	    }
	}
    }

    method GetOffset {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Offset
    }

    method ComputeNormalizedScaleFactor {} {
	# Look at both min and max in dataset for scaling since
	# large negative values can will get scaled below plane

	set plane_min [lindex $TableRange 0]
	set plane_max [lindex $TableRange 1]

	set data_range [expr abs($plane_max - $plane_min)]

	# There are problems with OpenGL operating near machine epsilon
	# so for small differences force a larger data range.
	set absmax [expr abs($plane_min)]
	if {$absmax < [expr abs($plane_max)]} {
	    set absmax [expr abs($plane_max)]
	}

	if {[expr $data_range / $absmax] < [expr 100.0*$::cvis::math::FLT_EPSILON]} {
	    set data_range [expr $absmax*100.0*$::cvis::math::FLT_EPSILON]
	}

	scan [$Reader GetBounds] "%f %f %f %f %f %f" \
		RealBounds(min,x) RealBounds(max,x) \
		RealBounds(min,y) RealBounds(max,y) \
		RealBounds(min,z) RealBounds(max,z)

	set max_bounds [expr -$::cvis::math::FLT_MAX]
	foreach i "x y z" {
	    set t [expr $RealBounds(max,$i) - $RealBounds(min,$i)]
	    if {$max_bounds < $t} { 
		set max_bounds $t
	    }
	}

	set NormalizedScaleFactor [expr ($max_bounds / $data_range)]
    }

    method ExecuteSetScaleFactor {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set offset [expr $Offset - [lindex $TableRange 0]]
	switch $SlicePlane {
	    "x" {
		set scale "$ScaleFactor 1 1"
		set actor_position "[expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0 0"
	    }
	    "y" {
		set scale "1 $ScaleFactor 1"
		set actor_position "0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)] 0"
	    }
	    "z" {
		set scale "1 1 $ScaleFactor"
		set actor_position "0 0 [expr ($offset*$ScaleFactor*$NormalizedScaleFactor)]"
	    }
	}

	if [string length $Reader] { 
	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    while {[string length $item]} {
		if { [string length [info commands $this.$item.planeActor]] } {
		    eval $this.$item.planeActor SetScaleFactor $scale
		    eval $this.$item.planeActor SetPosition $actor_position
		}
		set item [$collection GetNextItem]
	    }

	    if {$BoundaryEdgesInterior} {
		set collection [$Reader GetCollection]
		$collection InitTraversal
		set patch 0
		set item [$collection GetNextItem]
		while {[string length $item]} {
		    if { [string length [info commands \
			eval $this.$item.boundaryActor SetPosition $actor_position
		    }
		    incr patch
		    set item [$collection GetNextItem]
		}
	    } {
		for {set level 0} \
			{$level < [[$Reader GetReader] GetNumberOfLevels]} \
			{incr level} {
		    if {$LevelHasData($level)} {
			if {[string length [info commands \
			    eval $this.$level.boundaryActor SetPosition $actor_position
			}
		    }
		}
	    }
	}
    }

    method SetScaleFactor {factor} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$ScaleFactor != $factor} {

	    set ScaleFactor $factor

	    ComputeNormalizedScaleFactor
	
	    ::cvis::queue::Add "$this ExecuteSetScaleFactor"
	}
    }

    
    method GetScaleFactor {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $ScaleFactor
    }

    method SetTableRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$TableRange != "$min $max"} {

	    set TableRange "$min $max"

	    ComputeNormalizedScaleFactor

	    ::cvis::queue::Add "$this ExecuteSetTableRange"
	}
    }

    method ExecuteSetTableRange {} {
	if {!$DataCreateMode} {
	    set collection [$Reader GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
		while {[string length $item]} {
		    set item [$collection GetNextItem]
		}
	    }
	    
	    ExecuteSetScaleFactor
	}
    }

    method GetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $TableRange
    }

    method GetScalarRange {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {[string length $Reader]} {
	    return [$Reader GetScalarRange]
	} {
	    return "0 1"
	}
    }

    method SetLookupTable {lut} {
	set tracervar [::itcl::local cvisTrace #auto]
	$lut AddDepends $this
	set LookupTable $lut
	UpdateLookupTable
    }

    method UpdateLookupTable {} {
	set tracervar [::itcl::local cvisTrace #auto]

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

    variable Depends ""

    method AddDepends { depend } {
	set tracervar [::itcl::local cvisTrace #auto]

	lappend Depends $depend
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]


	if {[string length $Reader]} {
	    set Bounds [$Reader GetBounds]
	}


	# SGS is this really needed?
	if {1} {
	    # Update user interface elements
	    if {[string length $Interface]} {
		$Interface Update
		# Get values that might have changed based on UI control
		set Slice [$Interface GetSliceInScreenCoord]
	    }
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	if {$DataExists} {
	    set DataExists 0
	    ::cvis::queue::Add "$this ExecuteDataDelete"
	}

	if {$BoundaryExists} {
	    set BoundaryExists 0
	    ::cvis::queue::Add "$this ExecuteBoundaryDelete"
	}

	# Reset the collection information
	::cvis::queue::Add "$this ExecuteSetPreviousCollection"

	if {$DataIsVisible} {
	    set InDataCreateMode 1
	    ::cvis::queue::Add "$this ExecuteDataCreate"
	}

	if {$BoundaryEdgesVisible} {
	    ::cvis::queue::Add "$this ExecuteBoundaryCreate"
	}

	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}
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
}


