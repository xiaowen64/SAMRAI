
class cvisDataBrowser {
    inherit cvisDisplayObject

    variable PreviousCollection ""
    variable PreviousLevel 0

    variable ScaleFactor 1.0
    variable ScalarRange "0.0 1.0"

    variable LookupTable ""

    variable Interface ""

    variable AxisScale "1.0 1.0 1.0"

    variable BoundaryEdgesVisible 1
    variable BoundaryEdgeWidth 3

    variable DataVisible 1
    variable LabelVisible 1
    variable RenderWindowVisible 1

    variable Input ""

    variable SlicePlane "z"
    variable PickPosition
    variable PickWindowSize

    variable PickWindow

    variable DegreesToRadians 0.0

    variable Bounds 
    variable Delta 
    variable Scaling
    
    variable GlobalScaling 1.0

    method SetScaling {x y z} {
	set Scaling(x) $x
	set Scaling(y) $y
	set Scaling(z) $z
    }

    method SetDelta {x y z} {
	set Delta(x) $x
	set Delta(y) $y
	set Delta(z) $z
    }

    method SetCellBoundaryVisibility {visible} {

	set BoundaryEdgesVisible $visible 

  	if {$RenderWindowVisible} {
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    while {[string length $item]} {
			$BoundaryEdgesVisible
		set item [$Input GetNextItem]
	    }
	    }
	}
    }

    variable LabelFormat "%1.2e"
    method SetLabelFormat {format} {
	set LabelFormat $format
  	if {$RenderWindowVisible} {
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    set patch 0
	    while {[string length $item]} {
		
		set item [$Input GetNextItem]
		incr patch
	    }
	    
	}
    }
    method GetLabelFormat {} {
	return $LabelFormat
    }

    method GetCellBoundaryVisibility { } {
	return $BoundaryEdgesVisible
    }

    method SetDataVisibility {visible} {

	set DataVisible $visible 

  	if {$RenderWindowVisible} {
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    while {[string length $item]} {
			$DataVisible
		set item [$Input GetNextItem]
	    }
	    }
	}
    }

    method GetDataVisibility { } {
	return $DataVisible
    }

    method SetLabelVisibility {visible} {

	set LabelVisible $visible 

	if {$RenderWindowVisible} {
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    while {[string length $item]} {
			$LabelVisible
		set item [$Input GetNextItem]
	    }
	    }
	}
    }

    method GetLabelVisibility { } {
	return $LabelVisible
    }

    method SetPickWindow { window } {
	set PickWindow $window
	
	    
	    
	    




	    set NewPicker [$PickWindow GetPicker]
	    $NewPicker SetEndPickMethod  "$this PickMethod"
	}
    }

    variable NewPicker ""

    
    method PickMethod { } {
	if { [$NewPicker GetCellId] < 0 } {
	    if {$RenderWindowVisible} {
	    }
	} {
	    
	    scan [$NewPicker GetPickPosition] "%f %f %f" \
		    pick_position(x) pick_position(y) pick_position(z)

	    scan [$NewPicker GetSelectionPoint] "%f %f" \
		    selectionpoint(x) selectionpoint(y)

	    if {$RenderWindowVisible} {
	    
		set camera [$renderer GetActiveCamera]
		scan [$camera GetFocalPoint] "%f %f %f" \
			focal_point(x) focal_point(y) focal_point(z)
		scan [$camera GetPosition] "%f %f %f" \
			position(x) position(y) position(z)
		
		# Based on the plane we are viewing set the camera
		# in the databrowser window
		switch $SlicePlane { 
		    "x" {
				$pick_position(z) $pick_position(y) 
			
			# Picker in DataBrowser window is in 2D space
				0 $pick_position(y) $pick_position(z)

		    }
		    
		    "y" {
				$pick_position(x) $pick_position(z) 
			
			# Picker in DataBrowser window is in 2D space
				$pick_position(x) 0 $pick_position(z) 
		    }
		    
		    "z" {
				$pick_position(x) $pick_position(y)
			
			# Picker in DataBrowser window is in 2D space
				$pick_position(x) $pick_position(y) 0

			#******************************
			#******************************
			#******************************
			set PickPosition(x) $pick_position(x)
			set PickPosition(y) $pick_position(y)
			#******************************
			#******************************
			#******************************
		    }
		}

		$camera ComputeViewPlaneNormal
		
		# Set the size of the picker
		set distance [$camera GetDistance] 
		set view_angle [$camera GetViewAngle]
		set scale [$camera GetParallelScale]
		set view_angle [expr $DegreesToRadians * $view_angle]
		
		set width [expr $scale * 0.05]
		

	    }
	    
	    # Selector in render window is in 3D space
	    set camera [[$PickWindow GetRenderer] GetActiveCamera]
	    
		    $pick_position(x) $pick_position(y) $pick_position(z)

	    # Set the size of the selector
	    set distance [$camera GetDistance] 
	    set view_angle [$camera GetViewAngle]
	    set view_angle [expr $DegreesToRadians * $view_angle]
	    set scale [$camera GetParallelScale]

	    set width [expr $distance * tan($view_angle / 2.0) * 0.01]

	    
	    # Text message for pick location
	    
		    $selectionpoint(x) $selectionpoint(y)

	}
    }
    
    method SetSlicePlane {plane} {
	if {$SlicePlane != $plane}  {
	    set SlicePlane $plane
	    if {$RenderWindowVisible} {
	    }
	}
    }

    method GetSlicePlane { } {
	return $SlicePlane
    }

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	vtkMath math 
	set DegreesToRadians [math DegreesToRadians]
	math Delete

	set PickPosition(x) 4.0
	set PickPosition(y) 8.8

	set PickWindowSize(x) 4
	set PickWindowSize(y) 4

	foreach plane "x y z" { 
	    set Delta($plane)      0.1
	    set Bounds(min,$plane) 0.0
	    set Bounds(max,$plane) 1.0
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

    method SetInput {input} {
	set Input $input
	Update
    }

    method SetScalarRange {min max} {
	set ScalarRange "$min $max"
    }

    method GetScalarRange {} {
	return $ScalarRange 
    }

    method SetLookupTable {lut} {
	set LookupTable $lut
    }

    variable Depends ""

    method AddDepends { depend } {
	lappend Depends $depend
    }

    method SetCellBoundaryEdgeWidth { width } {
	set BoundaryEdgeWidth $width

  	if {[string length $Input]} {
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    while {[string length $item]} {
		set item [$Input GetNextItem]
	    }
	    }

	}
    }

    method GetCellBoundaryEdgeWidth { } {
	return $BoundaryEdgeWidth
    }

    method GetBounds { } {
	ComputeBounds
	return "$Bounds(min,x) $Bounds(max,x) $Bounds(min,y) $Bounds(max,y) $Bounds(min,z) $Bounds(max,z)"
    }

    method ComputeBounds { } {
  	if {[string length $Input]} {
	    
	    foreach plane "x y z" { 
		set Bounds(min,$plane) $::cvis::math::FLT_MAX
		set Bounds(max,$plane) [expr -$::cvis::math::FLT_MAX]
	    }
	    
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    set patch 0
	    while {[string length $item]} {
		scan [$item GetBounds] "%f %f %f %f %f %f" \
			new_real(min,x) new_real(max,x) \
			new_real(min,y) new_real(max,y) \
			new_real(min,z) new_real(max,z)
		
		foreach plane "x y z" {
		    if { $Bounds(min,$plane) > $new_real(min,$plane) } {
			set Bounds(min,$plane) $new_real(min,$plane)
		    }
		
		    if { $Bounds(max,$plane) < $new_real(max,$plane) } {
			set Bounds(max,$plane) $new_real(max,$plane)
		    }
		}
		
		set item [$Input GetNextItem]
	    }
	    
	} 
    }

    method GetSpacing { } {
	ComputeSpacing
	return "$Delta(x) $Delta(y) $Delta(z)"
    }

    method ComputeSpacing { } {
  	if {[string length $Input]} {
	    
	    foreach plane "x y z" { 
		set Delta($plane) $::cvis::math::FLT_MAX
	    }
	    
	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    set patch 0
	    while {[string length $item]} {
		scan [$item GetSpacing] "%f %f %f" \
			new_real(x) \
			new_real(y) \
			new_real(z)
		
		foreach plane "x y z" {
		    if { $Delta($plane) > $new_real($plane) } {
			set Delta($plane) $new_real($plane)
		    }
		}
		
		set item [$Input GetNextItem]
	    }
	    
	}
    }

    method Update {} {
	global NumBorderColors
	global BorderColors
	# Get rid of the old items
	foreach item $PreviousCollection {




	    }
	}
	
	if { $RenderWindowVisible } {

	    # If the data browser rendering window does not exist 
	    # create it



			SetParallelProjection 1
	    }

	    #****************************************************
	    #****************************************************
	    #****************************************************

	    ComputeBounds
	    ComputeSpacing


	    }

	    set NumCells 0
	    
	    #***************************************************
	    #***************************************************
	    #***************************************************

	    $Input InitTraversal
	    set item [$Input GetNextItem]
	    set patch 0
	    while {[string length $item]} {

		puts "$item *****************************************"
		lappend PreviousCollection $item

		$item ComputeBounds
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

		#***************************************************
		#***************************************************
		#***************************************************
		switch $SlicePlane {
		    "z" {

			# Compute the index into the colums/rows from
			# the real values
			set origin_index(x) \
				[expr int(($origin(x) - $Bounds(min,x)) / \
				$Delta(x))]

			puts "$item Origin $origin(x) $origin(y)"

			puts "$item Bounds $real(min,x) $real(max,x) $real(min,y) $real(max,y)"

			set origin_index(y) \
				[expr int(($origin(y) - $Bounds(min,y)) / \
				$Delta(x))]

			puts "$item OriginIndex $origin_index(x) $origin_index(y)"

			set end_index(x) \
				[expr int(($real(max,x) - $Bounds(min,x)) / \
				$Delta(x))]

			set end_index(y) \
				[expr int(($real(max,y) - $Bounds(min,y)) / \
				$Delta(y))]

			puts "$item end $end_index(x) $end_index(y)"

			set step_index(x) \
				[expr int( $spacing(x) /$Delta(x))]

			set step_index(y) \
				[expr int( $spacing(y) /$Delta(y))]

			puts "$item step $step_index(x) $step_index(y)"

			set scalars [[$item GetCellData] GetScalars]

			set skip 0

			puts "[expr $PickPosition(x) + $PickWindowSize(x) * $Delta(x) ]"
			puts "[expr $PickPosition(x) - $PickWindowSize(x) * $Delta(x) ]"
			if {[expr $PickPosition(x) + \
				$PickWindowSize(x) * $Delta(x) ] \
				< $real(min,x)} {
			    set skip 1
			}

			if {[expr $PickPosition(x) - \
				$PickWindowSize(x)] * $Delta(x) \
				> $real(max,x)} {
			    set skip 1
			}

			if {[expr $PickPosition(y) + \
				$PickWindowSize(y) * $Delta(y) ] \
				< $real(min,y)} {
			    set skip 1
			}

			if {[expr $PickPosition(y) - \
				$PickWindowSize(y)] * $Delta(y) \
				> $real(max,y)} {
			    set skip 1
			}


			if {!$skip} {

			    puts "$item Drawing for ($PickPosition(x), $PickPosition(y))"

			    set pickwindowbounds [expr round($PickPosition(y) / $Delta(y)) - $PickWindowSize(y)]
			    puts "WindowBounds min,y $pickwindowbounds"
			    if {$pickwindowbounds < $origin_index(y)} {
				set start_j 0
			    } {
				set start_j [expr $pickwindowbounds - $origin_index(y)] 
				# Convert to local index space
				set start_j [expr int($start_j * $Delta(y)/$spacing(y))]
			    }

			    set pickwindowbounds [expr int($PickPosition(y) / $Delta(y)) + $PickWindowSize(y)]
			    puts "WindowBounds max,y $pickwindowbounds"
			    puts "$pickwindowbounds -  $origin_index(y) + ($dim(y) - 1)"
			    if {$pickwindowbounds > [expr $origin_index(y) + ($dim(y) - 1)*$Delta(y)/$spacing(y)] } {
				set end_j [expr $dim(y) - 1]
			    } {
				set end_j [expr $pickwindowbounds - $origin_index(y)]
				set start_j [expr int($start_j * $Delta(y)/$spacing(y))]
			    }
			    
			    puts "Extracting Region $start_j to $end_j"



			    set pickwindowbounds [expr round($PickPosition(x) / $Delta(x)) - $PickWindowSize(x)]
			    puts "WindowBounds min,x $pickwindowbounds"
			    if {$pickwindowbounds < $origin_index(x)} {
				set start_i 0
			    } {
				set start_i [expr $pickwindowbounds - $origin_index(x)] 
				# Convert to local index space
				set start_i [expr int($start_i * $Delta(x)/$spacing(x))]
			    }

			    set pickwindowbounds [expr int($PickPosition(x) / $Delta(x)) + $PickWindowSize(x)]
			    puts "WindowBounds max,x $pickwindowbounds"
			    puts "$pickwindowbounds -  $origin_index(x) + ($dim(x) - 1)"
			    if {$pickwindowbounds > [expr $origin_index(x) + ($dim(x) - 1)*$Delta(x)/$spacing(x)] } {
				set end_i [expr $dim(x) - 1]
			    } {
				set end_i [expr $pickwindowbounds - $origin_index(x)]
				set start_i [expr int($start_i * $Delta(x)/$spacing(x))]
			    }

			    puts "Extracting Region $start_i to $end_i"
			    
			    #		    for {set j 0} {$j < [expr $dim(y) - 1]} {incr j} 
			    # for {set i 0} {$i < [expr $dim(x) - 1]} {incr i} 

			    for {set j $start_j} {$j < $end_j} {incr j} {
				for {set i $start_i} {$i < $end_i} {incr i} {
				    set index_x [expr $origin_index(x) + $i * $step_index(x)]
				    set index_y [expr $origin_index(y) + $j * $step_index(y)]
				    
				    set scalar_offset [expr int($j * ($dim(x)-1) + $i)]
				    
				    frame $frame.cellframe$NumCells -borderwidth 1 -relief ridge -width 40 -height 40
				    label $frame.cellframe$NumCells.label \
					    -text "[format $LabelFormat [$scalars GetScalar $scalar_offset]]"
				    pack $frame.cellframe$NumCells.label

				    grid $frame.cellframe$NumCells \
					    -row [expr int($Bounds(max,y)/ $Delta(y) - $index_y - ($step_index(y) -1))] \
					    -rowspan $step_index(y) \
					    -column $index_x -columnspan  $step_index(x) \
					    -sticky news
				    
				    incr NumCells
				}
			    }
			}
		    }
		}
		#***************************************************
		#***************************************************
		#***************************************************


		# SGS This should be based on the slice
		switch $SlicePlane { 
		    "x" {
			$item SetOrigin 0 $origin(y) $origin(z) 
		    }
		    "y" {
			$item SetOrigin $origin(x) 0 $origin(z)
		    }
		    "z" {
			$item SetOrigin $origin(x) $origin(y) 0
		    }
		}

		scan [$item GetOrigin] "%f %f %f" \
			origin(x) origin(y)  origin(z)
		
		vtkStructuredPointsGeometryFilter \



		##---------------------------------------------------------
		# Cell Boundaries 
		##---------------------------------------------------------

		


		# SGS how do we get the level information ?
		set level 0
			$BorderColors([expr $level % $NumBorderColors])
		
			$BoundaryEdgesVisible
		
		

		
		##---------------------------------------------------------
		# Convert to Cell to Node centered data
		# for labeling
		##---------------------------------------------------------

		# SGS make this a programmable filter!
		    
		#SGS If dim becomes 1 do we need to do something
		# special?  How do we intersect?
		# Shrink size
		# if one D can bypass the cut and display
		# if it is in the same plane otherwise we get a line
		foreach plane "x y z" {
		    if { $plane != $SlicePlane } {
			incr dim($plane) -1
			set origin($plane) \
				[expr $origin($plane) + $spacing($plane)/2.0]
		    }
		}



			$spacing(x) $spacing(y) $spacing(z)
			$origin(x) $origin(y)  $origin(z)
		    
			[[$item GetCellData] GetScalars]

		DataBrowser.$item.browserData Update
		# SGS this is very hacky 
			SetScalars [[$item GetCellData] GetScalars]
		#---------------------------------------------------------
		# Label data
		#---------------------------------------------------------

		    



		#---------------------------------------------------------
		# Color the Cells by value
		#---------------------------------------------------------
		
			SetImmediateModeRendering $::cvis::options::LOWMEM

			$ScalarRange
		
		if {[string length $LookupTable] } { 
			    [$LookupTable GetLookupTable]
		}
		

		

		set item [$Input GetNextItem]
		incr patch
	    }

	    
	    }
	    # 
	    # SGS What is a good scaling value?

	    


	    ComputeBounds
	    ComputeSpacing

	    
	    switch $SlicePlane { 
		"x" {
		    # Z x Y Plane
			    $Bounds(max,y) $Bounds(min,y) $Delta(y) "Y"

			    $Bounds(min,z) $Bounds(max,z) $Delta(x) "Z"

		    set base_zoom [expr $Bounds(max,z) - $Bounds(min,z)]
	    
		    if {$base_zoom < [expr $Bounds(max,y) - $Bounds(min,y)]} {
			set base_zoom [expr $Bounds(max,y) - $Bounds(min,y)]
		    }
	    
			    1.0 100.0 1.0  $base_zoom \
			    "Zoom X"

			    [expr ($Bounds(max,z) - $Bounds(min,z)) / 2.0 ] \
			    [expr ($Bounds(max,y) - $Bounds(min,y)) / 2.0 ]

		}
		"y" {
		    # X x Z Plane
			    $Bounds(max,z) $Bounds(min,z) $Delta(y) "Z"

			    $Bounds(min,x) $Bounds(max,x) $Delta(x) "X"

		    set base_zoom [expr $Bounds(max,x) - $Bounds(min,x)]
	    
		    if {$base_zoom < [expr $Bounds(max,z) - $Bounds(min,z)]} {
			set base_zoom [expr $Bounds(max,z) - $Bounds(min,z)]
		    }
	    
			    1.0 100.0 1.0  $base_zoom \
			    "Zoom Y"

			    [expr ($Bounds(max,x) - $Bounds(min,x)) / 2.0 ] \
			    [expr ($Bounds(max,z) - $Bounds(min,z)) / 2.0 ]

		}
		"z" {
			    $Bounds(max,y) $Bounds(min,y) $Delta(y) "Y"

			    $Bounds(min,x) $Bounds(max,x) $Delta(x) "X"

		    set base_zoom [expr $Bounds(max,x) - $Bounds(min,x)]
	    
		    if {$base_zoom < [expr $Bounds(max,y) - $Bounds(min,y)]} {
			set base_zoom [expr $Bounds(max,y) - $Bounds(min,y)]
		    }
	    
			    1.0 100.0 1.0  $base_zoom \
			    "Zoom Z"

			    [expr ($Bounds(max,x) - $Bounds(min,x)) / 2.0 ] \
			    [expr ($Bounds(max,y) - $Bounds(min,y)) / 2.0 ]

		}
	    }

	} {
	    # Rendering window exists get rid of it

		}

	    }
	}
    }
}
