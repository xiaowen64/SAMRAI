##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisBoundingSlice.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Bounding Box
##

class cvisBoundingSlice {
    inherit cvisDisplayObject
    
    variable Reader ""
    
    variable PreviousLevel 0
    
    variable PatchEdgesVisible 1
    variable PatchEdgeWidth 1.0
    
    variable Interface ""
    
    variable LevelHasData 

    variable Level

    variable NumPatchColors 8
    variable PatchColors
    variable PatchVisible

    variable NumPatches 0

    variable PatchExists 0

    variable Slice 0.0
    variable SlicePlane "z"

    method GetSlice {} {
	return $Slice
    }

    method SetSlice {value} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$Slice == $value} {
	    return
	}

	set Slice $value
	ExecuteSetSlice
    }

    method ExecuteSetSlice {} {
	set tracervar [::itcl::local cvisTrace #auto]
	# Schedule CPU intensive tasks that accomplish the update
	if {$PatchExists} {
	    set PatchExists 0
	    ::cvis::queue::Add "$this ExecutePatchDelete"
	}

	# SGS Visiblity control via visibility or creation/destruction?
	::cvis::queue::Add "$this ExecutePatchCreate"
    }

    method GetSlicePlane {} {
	return $SlicePlane
    }
    method SetSlicePlane {plane} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {$SlicePlane == $plane} {
	    return
	}

	set SlicePlane $plane
	ExecuteSetSlice
    }

    method SetPatchColor {level r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$PatchColors($level) != "$r $g $b"} {
	    set PatchColors($level) "$r $g $b"
	    if {[string length $Reader]} {
		set patch 0
		while {$patch < $NumPatches} {
		    if { $Level($patch) == $level } {
			if {[string length \
				[info commands $this.$patch.patchActor]]} {
			    eval [$this.$patch.patchActor GetProperty] \
				    SetColor $PatchColors([expr $Level($patch) \
				    % $NumPatchColors])
			}
		    }
		    incr patch
		}
	    }
	}
    }

    method SetPatchVisibility {level visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$PatchVisible($level) != $visible} {
	    set PatchVisible($level) $visible
	    if {[string length $Reader]} {
		set patch 0
		while {$patch < $NumPatches} {
		    if { $Level($patch) == $level } {
			if {[string length \
				[info commands $this.$patch.patchActor]]} {
			    $this.$patch.patchActor SetVisibility \
				    [expr $PatchVisible([expr $Level($patch) \
				    % $NumPatchColors]) & \
				    $PatchEdgesVisible]
			}
		    }
		    incr patch
		}
	    }
	}
    }

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	set NumPatchColors 8
	# green
	set PatchColors(0) "0 1 0"
	# yellow
	set PatchColors(1) "1 1 0"
	# light blue
	set PatchColors(2) "0 1 1"
	# purple
	set PatchColors(3) "1 0 1"
	# pink 
	set PatchColors(4) "1 0.75 0.80"
	# orange
	set PatchColors(5) "1 0.65 0"
	# red
	set PatchColors(6) "1 0 0"
	# blue
	set PatchColors(7) "0 0 1"
	
	set PatchVisible(0) "1"
	set PatchVisible(1) "1"
	set PatchVisible(2) "1"
	set PatchVisible(3) "1"
	set PatchVisible(4) "1"
	set PatchVisible(5) "1"
	set PatchVisible(6) "1"
	set PatchVisible(7) "1"
    }

    method GetVisibility {} {
	return $PatchEdgesVisible
    }
    
    method SetVisibility {vis } {
	set tracervar [::itcl::local cvisTrace #auto]
	if { $PatchEdgesVisible != $vis } {
	    set PatchEdgesVisible $vis
	    ExecuteSetVisibility
	}
    }
    
    method ExecuteSetVisibility { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    set patch 0
	    while {$patch < $NumPatches} {
		set level [$Reader GetPatchLevel $patch]
		if {[string length \
			[info commands $this.$patch.patchActor]]} {
		    $this.$patch.patchActor SetVisibility \
			    [expr $PatchVisible([expr $Level($patch) \
				    % $NumPatchColors]) & \
			    $PatchEdgesVisible]
		}
		incr patch
	    }
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

    method GetPatchEdgeWidth {} {
	return $width
    }
    
    method SetPatchEdgeWidth { width } {
	set tracervar [::itcl::local cvisTrace #auto]

	set PatchEdgeWidth $width
	
	if {[string length $Reader]} {
	    set patch 0
	    while {$patch < $NumPatches} {
		if {[string length \
			[info commands $this.$patch.patchActor]]} {
		    [$this.$patch.patchActor GetProperty] SetLineWidth \
			    $PatchEdgeWidth
		}
		incr patch
	    }
	}
    }
    
    method SetInput { reader } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Reader $reader
	$reader AddOutput $this
    }
    
    # Need since we are Output of file reader
    method SetScalarRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]
    }
    
    method SetInterface { interface } {
	set tracervar [::itcl::local cvisTrace #auto]
	set Interface $interface
	$interface Update
    }

    method ExecutePatchDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]
	set PatchExists 0
	set patch 0
	while {$patch < $NumPatches} {
	    if {[string length \
		    [info commands $this.$patch.outlineMapper]]} {
		$RenderWindow RemoveActor $this.$patch.patchActor
		::itcl::delete object $this.$patch.patchActor 
		$this.$patch.outlineMapper Delete
		$this.$patch.outlineSource Delete
	    }
	    incr patch
	}
    }

    method ExecutePatchCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set PatchExists 1
	if {[string length $Reader]} {

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

	    # Convert from real to screen coordinates
	    set scaling [$Reader GetScaling]

	    set NumPatches [$Reader GetNumberOfPatchBoundaries]
	    
	    set patch 0
	    while {$patch < $NumPatches} {

		set line [[$Reader GetPatchBoundary $patch] GetBounds]

		scan $line "%f %f %f %f %f %f" \
			bounds(min,x) bounds(max,x) \
			bounds(min,y) bounds(max,y) \
			bounds(min,z) bounds(max,z)

		if { $Slice >= $bounds(min,$SlicePlane) && \
			$Slice <= $bounds(max,$SlicePlane) } {

		    set level [$Reader GetPatchBoundaryLevel $patch]
		
		    # Only draw in the plane
		    set bounds(min,$SlicePlane) [expr $Slice - 10*$jitter]
		    set bounds(max,$SlicePlane) [expr $Slice + 10*$jitter]

		    set Level($patch) $level

		    vtkOutlineSource $this.$patch.outlineSource
		    eval $this.$patch.outlineSource SetBounds \
			    $bounds(min,x) $bounds(max,x) \
			    $bounds(min,y) $bounds(max,y) \
			    $bounds(min,z) $bounds(max,z)

		    vtkPolyDataMapper $this.$patch.outlineMapper
		    $this.$patch.outlineMapper SetImmediateModeRendering $::cvis::options::LOWMEM
		    $this.$patch.outlineMapper SetInput \
			    [$this.$patch.outlineSource GetOutput]
		    
		    cvisActor $this.$patch.patchActor
		    $this.$patch.patchActor SetMapper \
			    $this.$patch.outlineMapper
		    
		    $this.$patch.patchActor SetVisibility \
			    [expr $PatchVisible([expr $Level($patch) \
				    % $NumPatchColors]) & \
			    $PatchEdgesVisible]
		    
		    eval [$this.$patch.patchActor GetProperty] SetColor \
			    $PatchColors([expr $level % $NumPatchColors])
		    [$this.$patch.patchActor GetProperty] SetAmbient  1
		    [$this.$patch.patchActor GetProperty] SetDiffuse  0
		    [$this.$patch.patchActor GetProperty] SetSpecular 0
		    [$this.$patch.patchActor GetProperty] SetLineWidth \
			    $PatchEdgeWidth
		    
		    $this.$patch.patchActor SetPickable 0
		    
		    $RenderWindow AddActor $this.$patch.patchActor 
		}
		incr patch
	    }
	}
    }
    
    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Interface]} {
	    $Interface Update
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	# SGS do this
	$RenderWindow SetScaling

	# Schedule CPU intensive tasks that accomplish the update
	if {$PatchExists} {
	    set PatchExists 0
	    ::cvis::queue::Add "$this ExecutePatchDelete"
	}

	# SGS Visiblity control via visibility or creation/destruction?
	::cvis::queue::Add "$this ExecutePatchCreate"
    }
}



