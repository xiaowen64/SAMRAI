##
## File:        cvisDataBrowser.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: "Real" Data Browser
##

class cvisDataBrowser {
    inherit cvisDisplayObject

    variable ScalarRange "0.0 1.0"

    variable LookupTable ""

    variable Interface ""

    variable DataVisible 1
    variable LabelVisible 1
    variable RenderWindowVisible 0

    variable Reader ""

    variable SlicePlane "z"
    variable PickPosition
    variable PickWindowSize


    variable PickWindow

    variable DegreesToRadians 0.0

    variable Bounds 
    variable Scaling
    
    variable GlobalScaling 1.0

    variable DataLabelFormat "%1.2e"
    variable AxisLabelFormat "%1.2e"

    variable NewPicker ""

    variable NumCells 0
    variable GridFrames 

    method GetBounds { } { 
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetBounds]
	} {
	    return {0.0 1.0 0.0 1.0 0.0 1.0}
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

    method GetIndexScaling { } { 
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [[$Reader GetReader] GetIndexScaling]
	} {
	    return {1.0 1.0 1.0}
	}
    }

    method GetScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [$Reader GetScaling]
	} {
	    return 1.0
	}
    }

    method Shift {shift_x shift_y} {
	switch $SlicePlane {
	    "x" { 
		set sx "y"
		set sy "z"
	    }
	    "y" {
		set sx "x"
		set sy "z"
	    }
	    "z" {
		set sx "x"
		set sy "y"
	    }
	}

	set PickPosition($sx) [expr $PickPosition($sx) + \
		$shift_x * $Delta_dr($sx)]
	set PickPosition($sy) [expr $PickPosition($sy) + \
		$shift_y * $Delta_dr($sy)]

	if {[string length $Interface]} {
	    $Interface SetPickPoint \
		    $PickPosition(x) $PickPosition(y) $PickPosition(z)
	}

	if {[string length $RenderWindow] } { 
	    $RenderWindow SetPickPoint \
		    $PickPosition(x) $PickPosition(y) $PickPosition(z)
	}


	Update
    }

    method SetPickPoint {x y z} {
	set PickPosition(x) $x
	set PickPosition(y) $y
	set PickPosition(z) $z

	Update
    }

    method SetWindowSize {x y z} {
	set PickWindowSize(x) $x
	set PickWindowSize(y) $y
	set PickWindowSize(z) $z
	
	CreateWindow
	Update
    }

    method SetDataLabelFormat {format} {
	set DataLabelFormat $format
  	if {$RenderWindowVisible} {
	    LabelData
	}
    }

    method GetDataLabelFormat {} {
	return $DataLabelFormat
    }

    method SetAxisLabelFormat {format} {
	set AxisLabelFormat $format
  	if {$RenderWindowVisible} {
	    LabelAxis
	}
    }

    method GetAxisLabelFormat {} {
	return $AxisLabelFormat
    }

    method SetSlicePlane {plane} {
	if {$SlicePlane != $plane}  {
	    set SlicePlane $plane
	    if {$RenderWindowVisible} {
		CreateWindow
		Update
	    }
	}
    }

    method GetSlicePlane { } {
	return $SlicePlane
    }

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	set PickPosition(x) 0.5
	set PickPosition(y) 0.5
	set PickPosition(z) 0.5

	set PickWindowSize(x) 4
	set PickWindowSize(y) 4
	set PickWindowSize(z) 4

	foreach plane "x y z" { 
	    set Bounds(min,$plane) 0.0
	    set Bounds(max,$plane) 1.0
	    set Delta_sc($plane) 0.1
	}
    }

    method SetRenderWindowVisibility {visibility} {
	if { $visibility != $RenderWindowVisible } {
	    set RenderWindowVisible $visibility
	    Update
	}
    }
    
    method GetRenderWindowVisibility { } {
	return $RenderWindowVisible
    }


    method SetInput { reader } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Reader $reader
	$reader AddOutput $this
    }

    # SGS Do we need this?
    method SetScalarRange {min max} {
	set ScalarRange "$min $max"
    }

    # SGS Do we need this?
    method GetScalarRange {} {
	return $ScalarRange 
    }

    # SGS Do we need this?
    method SetLookupTable {lut} {
	set LookupTable $lut
    }

    variable Depends ""

    method AddDepends { depend } {
	lappend Depends $depend
    }

    method SetInterface { interface } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Interface $interface
	$interface Update
    }

    method ComputeIntersections {} {
	scan [$Reader GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z) 
	
	set Scaling [$Reader GetScaling]

	foreach plane {x y z} {
	    set Delta_sc($plane) [expr $Scaling * $IndexScaling($plane)]
	}
	
	switch $SlicePlane {
	    "x" { 
		set sx "y"
		set sy "z"
	    }
	    "y" {
		set sx "x"
		set sy "z"
	    }
	    "z" {
		set sx "x"
		set sy "y"
	    }
	}

	set numlevels [[$Reader GetReader] GetNumberOfRefinementLevels]

	set level [expr $numlevels - 1]
	while {$level > -1} {
	    foreach plane {x y z} {
		set Multiple($level,$plane) 1
	    }
	    incr level -1
	}

	set level [expr $numlevels - 1]
	while {$level > -1 } {
	    scan [[$Reader GetReader] GetRefinement $level] "%f %f %f" \
		    refinement(x) \
		    refinement(y) \
		    refinement(z) 
	    
	    set l [expr $level -1 ]
	    while {$l > -1} {
		foreach plane {x y z} {
		    set Multiple($l,$plane) \
			    [expr round($Multiple($l,$plane) * $refinement($plane))]
		}
		incr l -1
	    }
	    incr level -1
	}

	set max_level 0
	set level [expr $numlevels - 1]
	set continue 1
	while { $continue && ($level > -1)} {

	    set intersected ""
	    set levels ""

	    foreach  plane "x y z" {
		set Delta_dr($plane) \
			[expr $Delta_sc($plane) * $Multiple($level,$plane)]
		set PW_bounds_dg(min,$plane) \
			[expr round($PickPosition($plane) / $Delta_dr($plane)) \
			- $PickWindowSize($plane)]
		set PW_bounds_dg(max,$plane) \
			[expr round($PickPosition($plane) / $Delta_dr($plane)) \
			+ $PickWindowSize($plane) + 1]
	    }

	    set collection [[$Reader GetReader] GetCollection]
	    $collection InitTraversal
	    set item [$collection GetNextItem]
	    set patch 0
	    while {[string length $item]} {
		scan [$item GetBounds] "%f %f %f %f %f %f" \
			bounds_sr(min,x) bounds_sr(max,x) \
			bounds_sr(min,y) bounds_sr(max,y) \
			bounds_sr(min,z) bounds_sr(max,z)

		if {$PickPosition($SlicePlane) >= $bounds_sr(min,$SlicePlane) \
			&& $PickPosition($SlicePlane) <= \
			$bounds_sr(max,$SlicePlane)} {

		    scan [$item GetDimensions] "%d %d %d" \
			    dim_lg(x) dim_lg(y) dim_lg(z)

		    scan [$item GetSpacing] "%f %f %f" \
			    spacing_sr(x) spacing_sr(y)  spacing_sr(z)
		
		    scan [$item GetOrigin] "%f %f %f" \
			    origin_sr(x) origin_sr(y)  origin_sr(z)

		    # Compute the index into the colums/rows from
		    # the screen real values
		    foreach plane {x y z} {
			set origin_dg($plane) \
				[expr round($origin_sr($plane) / \
				$Delta_dr($plane))]
			set end_dg($plane) \
				[expr round($bounds_sr(max,$plane) / \
				$Delta_dr($plane))]
		    
			set step_dg($plane) \
			    [expr round( $spacing_sr($plane) / \
			    $Delta_dr($plane)) ]

			set pick_position_dg($plane) \
				[expr round($PickPosition($plane) / \
				$Delta_dr($plane))]

			set pick_position_lg($plane) \
			    [expr round(($PickPosition($plane) - \
			    $origin_sr($plane)) / $spacing_sr($plane))]
		    }

		    set skip 0
		
		    if {[expr $PickPosition($sx) + \
			    $PickWindowSize($sx) * $Delta_dr($sx) ] \
			    < $bounds_sr(min,$sx)} {
			set skip 1
		    }
		    
		    if {[expr $PickPosition($sx) - \
			    $PickWindowSize($sx) * $Delta_dr($sx) ] \
			    > $bounds_sr(max,$sx)} {
			set skip 1
		    }
		    
		    if {[expr $PickPosition($sy) + \
			    $PickWindowSize($sy) * $Delta_dr($sy) ] \
			    < $bounds_sr(min,$sy)} {
			set skip 1
		    }
		    
		    if {[expr $PickPosition($sy) - \
			    $PickWindowSize($sy) * $Delta_dr($sy) ] \
			    > $bounds_sr(max,$sy)} {
			set skip 1
		    }

		    if {$skip} {
		    } {
			lappend intersected $item
			lappend levels [[$Reader GetReader] \
				GetPatchLevel $patch]
		    }
		}

		set item [$collection GetNextItem]
		incr patch
	    }

	    if {[llength $intersected]}  {

		set max_level -1
		foreach item $intersected l $levels {
		    if {$max_level < $l} {
			set max_level $l
		    }
		}

		if {$max_level == $level} {
		    set continue 0
		} elseif { $max_level > $level } {
		    set continue 0
		    set max_level [expr $level +1]
		    set intersected $prev_intersected
		    set levels $prev_levels
		}
	    }
	    set prev_intersected $intersected
	    set prev_levels $levels
	    incr level -1	    
	}

	set Intersected $intersected
	set Level $max_level
	foreach  plane "x y z" {
	    set Delta_dr($plane) \
		    [expr $Delta_sc($plane) * $Multiple($Level,$plane)]
	    set PW_bounds_dg(min,$plane) \
		    [expr round($PickPosition($plane) / $Delta_dr($plane)) \
		    - $PickWindowSize($plane)]
	    set PW_bounds_dg(max,$plane) \
		    [expr round($PickPosition($plane) / $Delta_dr($plane)) \
		    + $PickWindowSize($plane) + 1]
	}

    }

    variable Intersected ""
    variable Level 0

    variable Multiple

    variable Delta_sc
    variable Delta_dr
    variable PW_bounds_dg

    method LabelAxis {} {

	    return
	}

	scan [$Reader GetIndexScaling] "%f %f %f" \
		IndexScaling(x) \
		IndexScaling(y) \
		IndexScaling(z) 
	
	set Scaling [$Reader GetScaling]
	
	switch $SlicePlane {
	    "x" { 
		set sx "y"
		set sy "z"
	    }
	    "y" {
		set sx "x"
		set sy "z"
	    }
	    "z" {
		set sx "x"
		set sy "y"
	    }
	}
	
	foreach  plane "$sx $sy" {
	    set PW_bounds_bg(min,$plane) [expr round($PickPosition($plane) / \
		    $Delta_dr($plane)) - $PickWindowSize($plane)]
	    set PW_bounds_bg(max,$plane) [expr round($PickPosition($plane) / \
		    $Delta_dr($plane)) + $PickWindowSize($plane) + 1]
	}

	#
	#
	# Label Plane
	$GridFrames(0,0).label configure -text \
		"[format "$SlicePlane=$AxisLabelFormat" [expr $PickPosition($SlicePlane) / $Scaling]]"

	
	
	$GridFrames(0,2).label configure -text "$sx" 
	for {set j $PW_bounds_bg(min,$sy)} { $j < $PW_bounds_bg(max,$sy)} \
		{incr j} {

	    set grid_row [expr $j - $PW_bounds_bg(min,$sy)]

	    # For display 0 0 is lower corner
	    set grid_row [expr 2*$PickWindowSize($sy)+2 - $grid_row]
	    
	    # make room for column labels
	    $GridFrames($grid_row,1).label configure -text \
		    "[format "$AxisLabelFormat" \
		    [expr ($j+0.5)*$Multiple($Level,$sy) \
		    * $IndexScaling($sy)]]"
	}

	$GridFrames(2,0).label configure -text "$sy" 
	for {set i $PW_bounds_bg(min,$sx)} \
		{ $i < $PW_bounds_bg(max,$sx)} {incr i} {
	    # reset to 0 origin
	    set grid_column [expr $i - $PW_bounds_bg(min,$sx)]
	    # make room for column labels
	    incr grid_column
	    incr grid_column
	    $GridFrames(1,$grid_column).label configure \
		    -text "[format "$AxisLabelFormat" \
		    [expr ($i+0.5)*$Multiple($Level,$sx) \
		    * $IndexScaling($sx)]]"
	}
    }

    method CreateWindow {} {
	}
	
	set NumAxis 0

	
	bind $frame <KeyPress-Left> "$this Shift -1 0"
	bind $frame <KeyPress-Right> "$this Shift 1 0"
	bind $frame <KeyPress-Up> "$this Shift 0 1"
	bind $frame <KeyPress-Down> "$this Shift 0 -1"
	
	switch $SlicePlane {
	    "x" { 
		set sx "y"
		set sy "z"
	    }
	    "y" {
		set sx "x"
		set sy "z"
	    }
	    "z" {
		set sx "x"
		set sy "y"
	    }
	}


	#
	# Label Plane
	#
	set GridFrames(0,0) $frame.axisframe_$NumAxis
	frame $GridFrames(0,0) \
		-borderwidth 1 -relief ridge 
		
	label $GridFrames(0,0).label -bg darkgrey 
	pack $GridFrames(0,0).label  -fill both -expand yes
	
	grid $GridFrames(0,0) -row 0 -column 0 -sticky news
	incr NumAxis

	#
	# Create labels for the screen X axis
	#
	
	set GridFrames(2,0) $frame.axisframe_$NumAxis
	frame $GridFrames(2,0) \
		-borderwidth 1 -relief ridge 
	label $GridFrames(2,0).label -text "$sx" -bg darkgrey 
	pack  $GridFrames(2,0).label  -fill both -expand yes
	grid $GridFrames(2,0) \
		-row 2 -column 0 \
		-rowspan [expr $PickWindowSize($sx) * 2 + 1] \
		-sticky news
	incr NumAxis
	
	for {set i 2} { $i < [expr $PickWindowSize($sx) * 2 + 3] } { incr i} {
	    set GridFrames(1,$i) $frame.axisframe_$NumAxis
	    frame $GridFrames(1,$i) \
		    -borderwidth 1 -relief ridge 
	    label $GridFrames(1,$i).label 
	    pack $GridFrames(1,$i).label  -fill both -expand yes
	    grid $GridFrames(1,$i) -row 1 -column $i -sticky news
	    incr NumAxis
	}
	
	#
	# Create labels for the screen Y axis
	#

	set GridFrames(0,2) $frame.axisframe_$NumAxis
	frame $GridFrames(0,2) \
		-borderwidth 1 -relief ridge 
	label $GridFrames(0,2).label -text "$sy" -bg darkgrey 
	pack  $GridFrames(0,2).label  -fill both -expand yes
	grid $GridFrames(0,2) \
		-row 0 -column 2 \
		-columnspan [expr $PickWindowSize($sy) * 2 + 1] \
		-sticky news
	incr NumAxis

	for {set j 2} { $j < [expr $PickWindowSize($sy) * 2 + 3] } { incr j} {
	    set GridFrames($j,1) $frame.axisframe_$NumAxis
	    frame $GridFrames($j,1) \
		    -borderwidth 1 -relief ridge 
	    label $GridFrames($j,1).label 
	    pack $GridFrames($j,1).label  -fill both -expand yes
	    grid $GridFrames($j,1) -row $j -column 1 -sticky news
	    incr NumAxis
	}

    }

    method LabelData {} {
	    return
	}
	
	
	if {$NumCells} {
	    incr NumCells -1
	    while {$NumCells} {
		destroy $frame.cellframe_$NumCells
		incr NumCells -1
	    }
	    destroy $frame.cellframe_$NumCells
	}
	
	set Scaling [$Reader GetScaling]
	
	switch $SlicePlane {
	    "x" { 
		set sx "y"
		set sy "z"
	    }
	    "y" {
		set sx "x"
		set sy "z"
	    }
	    "z" {
		set sx "x"
		set sy "y"
	    }
	}
	
	foreach item $Intersected {
	    scan [$item GetBounds] "%f %f %f %f %f %f" \
		    bounds_sr(min,x) bounds_sr(max,x) \
		    bounds_sr(min,y) bounds_sr(max,y) \
		    bounds_sr(min,z) bounds_sr(max,z)
	    
	    scan [$item GetDimensions] "%d %d %d" \
		    dim_lg(x) dim_lg(y) dim_lg(z)
		
	    scan [$item GetSpacing] "%f %f %f" \
		    spacing_sr(x) spacing_sr(y)  spacing_sr(z)
		
	    scan [$item GetOrigin] "%f %f %f" \
		    origin_sr(x) origin_sr(y)  origin_sr(z)
		
	    # Compute the index into the colums/rows from
	    # the real values

	    foreach plane {x y z} {
		set origin_dg($plane) \
			[expr round($origin_sr($plane) / $Delta_dr($plane))]

		set end_dg($plane) \
			[expr round($bounds_sr(max,$plane)/ $Delta_dr($plane))]

		set step_dg($plane) \
			[expr round( ($spacing_sr($plane) /$Delta_dr($plane)) )]
		    
		set pick_position_dg($plane) \
			[expr round($PickPosition($plane)/ $Delta_dr($plane))]
		    
		set pick_position_lg($plane) \
			[expr round(($PickPosition($plane) - \
			$origin_sr($plane)) / $spacing_sr($plane))]

	    }

	    set scalars [[$item GetCellData] GetScalars]
		
	    if {$PW_bounds_dg(min,$sx) < $origin_dg($sx)} {
		set start_i 0
	    } {
		set start_i [expr $PW_bounds_dg(min,$sx) - \
			$origin_dg($sx)] 
		# Convert to local index space
		set start_i [expr round($start_i * \
			$Delta_dr($sx)/$spacing_sr($sx))]
	    }
		    
	    if {$PW_bounds_dg(max,$sx) > $end_dg($sx)} {
		set end_i [expr $dim_lg($sx) - 1]
	    } {
		set end_i [expr $PW_bounds_dg(max,$sx) - \
			$origin_dg($sx)]
		set end_i [expr ceil($end_i * \
			$Delta_dr($sx)/$spacing_sr($sx))]
	    }
		    
	    if {$PW_bounds_dg(min,$sy) < $origin_dg($sy)} {
		set start_j 0
	    } {
		set start_j [expr $PW_bounds_dg(min,$sy) - \
				$origin_dg($sy)] 
		# Convert to local index space
		set start_j [expr round($start_j * \
			$Delta_dr($sy)/$spacing_sr($sy))]
	    }
		    
	    if {$PW_bounds_dg(max,$sy) > $end_dg($sy) } {
		set end_j [expr $dim_lg($sy) - 1]
	    } {
		set end_j [expr $PW_bounds_dg(max,$sy) - \
				$origin_dg($sy)]
		set end_j [expr int(ceil($end_j * \
				$Delta_dr($sy)/$spacing_sr($sy)))]
	    }

	   
	    for {set j $start_j} {$j < $end_j} {incr j} {
		for {set i $start_i} {$i < $end_i} {incr i} {
		    set x_bg [expr $origin_dg($sx) + \
			    $i * $step_dg($sx)]

		    set y_bg [expr $origin_dg($sy) + \
			    $j * $step_dg($sy)]
			    
		    switch $SlicePlane {
			"x" {
			    set scalar_offset [expr round($i * ($dim_lg(y)-1) * ($dim_lg(x)-1) + \
				    $j * ($dim_lg(x)-1) + $pick_position_lg($SlicePlane))]
			}
				
			"y" {
			    set scalar_offset [expr round($j * ($dim_lg(y)-1) * ($dim_lg(x)-1) + \
				    $pick_position_lg($SlicePlane) * ($dim_lg(x)-1) + $i)]
			}
				
			"z" {
			    set scalar_offset [expr round($pick_position_lg($SlicePlane) * ($dim_lg(y)-1) \
				    * ($dim_lg(x)-1) + \
				    $j * ($dim_lg(x)-1) +$i)]

			}
			
		    }
			    
			    
		    set grid_col $x_bg
		    set grid_col [expr $grid_col - $PW_bounds_dg(min,$sx)]
			    
		    set grid_row $y_bg			    
		    set grid_row [expr $PW_bounds_dg(max,$sy) - $grid_row]
		    set grid_row [expr $grid_row - ($step_dg($sy))]

		    # grid row 0-8, grid col 0-8
		    set step($sx) $step_dg($sx)
		    set step($sy) $step_dg($sy)
		    
		    if {$grid_row < 0} {
			# grid_col is negative, the part hanging off the bottom
			set step($sy) [expr $step_dg($sy) + $grid_row]
			set grid_row 0
		    }

		    if {$grid_col < 0} {
			# grid_col is negative, the part hanging off the side
			set step($sx) [expr $step_dg($sx) + $grid_col]
			set grid_col 0
		    }

		    # make room for the row/column labels
		    incr grid_row 2
		    incr grid_col 2 
			    
		    if { ($pick_position_dg($sx) >= $x_bg &&  \
			    $pick_position_dg($sx) <= \
			    [expr $x_bg + $step_dg($sx)-1]) && \
			    ($pick_position_dg($sy) >= $y_bg &&  \
			    $pick_position_dg($sy) <= \
			    [expr $y_bg + $step_dg($sy)-1]) } {
			set color "grey"
		    } {
			set color "white"
		    }
			    
		    if {$step($sy) > 0 && $step($sx) > 0} {
			set cellframe $frame.cellframe_$NumCells
			frame $cellframe -borderwidth 1 \
				-relief ridge 
			label $cellframe.label \
				-bg $color \
				-text "[format $DataLabelFormat [$scalars GetValue $scalar_offset]]"
			pack $cellframe.label -fill both -expand yes
			grid $cellframe \
				-row $grid_row  \
				-rowspan  $step($sy) \
				-column $grid_col  \
				-columnspan  $step($sx) \
				-sticky news
			incr NumCells
		    }
		}
	    }
	}
    }

    method ComputeBoundingBox {} {
	set small_delta [expr 10*$::cvis::math::FLT_EPSILON]
	switch $SlicePlane {
	    "x" {
		set bounds \
			"[expr $PickPosition(x)-$small_delta] [expr $PickPosition(x)+$small_delta]  [expr $PickPosition(y) - $Delta_dr(y)*($PickWindowSize(y)+0.5)] [expr $PickPosition(y) + $Delta_dr(y)*($PickWindowSize(y)+0.5)] [expr $PickPosition(z) - $Delta_dr(z)*($PickWindowSize(z)+0.5)] [expr $PickPosition(z) + $Delta_dr(z)*($PickWindowSize(z)+0.5)]"
	    }

	    "y" {
		set bounds \
			"[expr $PickPosition(x) - $Delta_dr(x)*($PickWindowSize(x)+0.5)] [expr $PickPosition(x) + $Delta_dr(x)*($PickWindowSize(x)+0.5)] [expr $PickPosition(y)-$small_delta] [expr $PickPosition(y)+$small_delta] [expr $PickPosition(z) - $Delta_dr(z)*($PickWindowSize(z)+0.5)] [expr $PickPosition(z) + $Delta_dr(z)*($PickWindowSize(z)+0.5)]"
	    }

	    "z" {
		set bounds \
			"[expr $PickPosition(x) - $Delta_dr(x)*($PickWindowSize(x)+0.5)] [expr $PickPosition(x) + $Delta_dr(x)*($PickWindowSize(x)+0.5)] [expr $PickPosition(y) - $Delta_dr(y)*($PickWindowSize(y)+0.5)] [expr $PickPosition(y) + $Delta_dr(y)*($PickWindowSize(y)+0.5)] [expr $PickPosition(z)-$small_delta] [expr $PickPosition(z)+$small_delta]"
	    }
	}

	eval $RenderWindow SetPickBounds $bounds
    }
    
    method Update {} {

	if {![string length $Reader]} {
	    return
	}

	# Update user interface elements
	if {[string length $Interface]} {
	    $Interface Update
	    # Get values that might have changed based on UI control
	    # set Slice [$Interface GetSliceInScreenCoord]
	}

	if { $RenderWindowVisible } {
		CreateWindow
	    }
	    ComputeIntersections
	    ComputeBoundingBox
	    LabelAxis
	    LabelData

	} { 
	    eval $RenderWindow SetPickBounds "0 0 0 0 0 0"
	    # Rendering window exists get rid of it
	    }
	}
    }
}

