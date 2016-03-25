##
## File:        cvisBoundingBox.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Bounding Box
##

class cvisBoundingBox {
    inherit cvisDisplayObject
    
    variable Reader ""
    
    variable PreviousLevel 0
    
    variable BoundaryEdgesVisible 1
    variable BoundaryEdgeWidth 1.0
    
    variable Interface ""
    
    variable LevelHasData 

    variable Level

    variable NumBorderColors 8
    variable BorderColors
    variable BorderVisible

    variable NumPatches 0

    variable BoundaryExists 0

    method SetBorderColor {level r g b} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BorderColors($level) != "$r $g $b"} {
	    set BorderColors($level) "$r $g $b"
	    if {[string length $Reader]} {
		set patch 0
		while {$patch < $NumPatches} {
		    if { $Level($patch) == $level } {
			eval [$this.$patch.boundaryActor GetProperty] \
				SetColor $BorderColors([expr $Level($patch) \
				% $NumBorderColors])
		    }
		    incr patch
		}
	    }
	}
    }

    method SetBorderVisibility {level visible} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$BorderVisible($level) != $visible} {
	    set BorderVisible($level) $visible
	    if {[string length $Reader]} {
		set patch 0
		while {$patch < $NumPatches} {
		    if { $Level($patch) == $level } {
			$this.$patch.boundaryActor SetVisibility \
				[expr $BorderVisible([expr $Level($patch) \
				% $NumBorderColors]) & \
				$BoundaryEdgesVisible]
		    }
		    incr patch
		}
	    }
	}
    }

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
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
    
    method SetVisibility {vis } {
	set tracervar [::itcl::local cvisTrace #auto]
	if { $BoundaryEdgesVisible != $vis } {
	    set BoundaryEdgesVisible $vis
	    ExecuteSetVisibility
	}
    }
    
    method ExecuteSetVisibility { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    set patch 0
	    while {$patch < $NumPatches} {
		set level [$Reader GetPatchLevel $patch]
		$this.$patch.boundaryActor SetVisibility \
			[expr $BorderVisible([expr $Level($patch) \
				% $NumBorderColors]) & \
			$BoundaryEdgesVisible]

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
    
    method SetBoundaryEdgeWidth { width } {
	set tracervar [::itcl::local cvisTrace #auto]
	set BoundaryEdgeWidth $width
	
	if {[string length $Reader]} {
	    set patch 0
	    while {$patch < $NumPatches} {
		[$this.$patch.boundaryActor GetProperty] SetLineWidth \
			$BoundaryEdgeWidth
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

    method ExecuteBoundaryDelete {} {
	set tracervar [::itcl::local cvisTrace #auto]
	set BoundaryExists 0
	set patch 0
	while {$patch < $NumPatches} {
	    if {[string length \
		$RenderWindow RemoveActor $this.$patch.boundaryActor
		::itcl::delete object $this.$patch.boundaryActor 
	    }
	    incr patch
	}
    }

    method ExecuteBoundaryCreate {} {
	set tracervar [::itcl::local cvisTrace #auto]

	set BoundaryExists 1
	if {[string length $Reader]} {

	    # Convert from real to screen coordinates
	    set scaling [$Reader GetScaling]

	    set NumPatches [$Reader GetNumberOfPatchBoundaries]
	    
	    set patch 0
	    while {$patch < $NumPatches} {

		set line [[$Reader GetPatchBoundary $patch] GetBounds]

		scan $line "%f %f %f %f %f %f" \
			bounds(min,x) bounds(min,y) bounds(min,z) \
			bounds(max,x) bounds(max,y) bounds(max,z) 

		set level [$Reader GetPatchBoundaryLevel $patch]
		
		foreach plane "x y z" {
		    set bounds(min,$plane) [expr $bounds(min,$plane) * $scaling]
		    set bounds(max,$plane) [expr $bounds(max,$plane) * $scaling]
		}
		
		set Level($patch) $level
		
			[[$Reader GetPatchBoundary $patch] GetOutput]
		
		cvisActor $this.$patch.boundaryActor
		$this.$patch.boundaryActor SetMapper \

		$this.$patch.boundaryActor SetVisibility \
			[expr $BorderVisible([expr $Level($patch) \
				% $NumBorderColors]) & \
			$BoundaryEdgesVisible]

		eval [$this.$patch.boundaryActor GetProperty] SetColor \
			$BorderColors([expr $level % $NumBorderColors])
		[$this.$patch.boundaryActor GetProperty] SetAmbient  1
		[$this.$patch.boundaryActor GetProperty] SetDiffuse  0
		[$this.$patch.boundaryActor GetProperty] SetSpecular 0
		[$this.$patch.boundaryActor GetProperty] SetLineWidth \
			$BoundaryEdgeWidth

		$this.$patch.boundaryActor SetPickable 0
		
		$RenderWindow AddActor $this.$patch.boundaryActor 
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
	if {$BoundaryExists} {
	    set BoundaryExists 0
	    ::cvis::queue::Add "$this ExecuteBoundaryDelete"
	}

	# SGS Visiblity control via visibility or creation/destruction?
	::cvis::queue::Add "$this ExecuteBoundaryCreate"
    }
}



