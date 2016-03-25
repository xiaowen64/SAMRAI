##
## File:        cvisExtractSlicePlane.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Slice plane
##

class cvisExtractSlicePlane {
    inherit cvisDisplayObject

    # Slice plane
    variable SlicePlane "z"
    variable Slice 0

    # Data slice
    variable DataExists 0

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

    variable OutPutCollection ""

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	vtkStructuredPointsCollection $this.collection
	set OutPutCollection $this.collection
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

    #=======================================================================
    # Data Slice
    #=======================================================================
    method ExecuteDataCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]
	
	if {$DataExists} {
	    return
	}

	set normal(x) {1 0 0}
	set normal(y) {0 1 0}
	set normal(z) {0 0 1}
	
	set DataExists 1

	$OutPutCollection RemoveAllItems
	
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

		    $spacing(x) $spacing(y) $spacing(z)

		    $origin(x) $origin(y) $origin(z)

		    [[$item GetCellData] GetScalars]


	    set item [$collection GetNextItem]
	}
	
	ExecuteDataSetSlice 
    }

    method ExecuteDataDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set DataExists 0
	foreach item $PreviousCollection {

	    if {[string length \
	    }
	}
    }

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
	# SGS is this really needed anymore for cell centered data?
	set jitter [expr ($real(max,$SlicePlane) - $real(min,$SlicePlane))\
		* $::cvis::math::BOUNDARY_JITTER ]
	
	if { $Slice < [expr $real(min,$SlicePlane) + $jitter] } {
	    set Slice [expr $Slice + $jitter]
	}
	
	if { $Slice > [expr $real(max,$SlicePlane) - $jitter] } {
	    set Slice [expr $Slice - $jitter]
	}

	$OutPutCollection RemoveAllItems
	
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


	    }
	    
	    set item [$collection GetNextItem]
	}
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

    method GetOutput {} {
	return $OutPutCollection
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
	
	if {$DataExists} {
	    set DataExists 0
	    ::cvis::queue::Add "$this ExecuteDataDelete"
	}
	
	# Reset the collection information
	::cvis::queue::Add "$this ExecuteSetPreviousCollection"
	::cvis::queue::Add "$this ExecuteDataCreate"
	
	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}
    }

    method SetDataVisibility {vis} {
    }

    method SetTextureMap {texturemap} {
    }

    method SetTextureMapInterpolate {texturemapinterpolate} {
    }

    method SetCellBoundaryWidth {width} {
    }

    method SetCellBoundaryVisibility {vis} {
    }

    method SetCellColor {level r g b} {
    }

    method SetCellVisibility {level visible} {
    }

    method SetPatchVisibility {level visible} {
    }

    method SetPatchBoundaryVisibility {vis} {
    }

    method SetPatchColor {level r g b} {
    }

    method SetPatchBoundaryWidth {width} {
    }

    method SetTableRange {min max} {
    }

    method SetAxisScale { scalex scaley scalez } {
    }

    method SetLookupTable {lut} {
    }

    method UpdateLookupTable {} {
    }
}


