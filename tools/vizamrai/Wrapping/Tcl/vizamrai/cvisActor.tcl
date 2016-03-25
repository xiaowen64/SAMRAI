##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisActor.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Actor
##

class cvisActorBase {
    variable Scale
    variable _vtkActor

    constructor {} {
    } {
	set Scale(x) 1.0
	set Scale(y) 1.0
	set Scale(z) 1.0
    }

    destructor {
	$_vtkActor Delete
    }

    method Update {} {

    }

    method SetScale {x y z} {
	set Scale(x) $x
	set Scale(y) $y
	set Scale(z) $z
	$_vtkActor SetScale $x $y $z
    }

    method SetMapper {mapper} {
	$_vtkActor SetMapper $mapper
    }

    method SetTexture {texture} {
	$_vtkActor SetTexture $texture
    }

    method SetVisibility {visible} {
	$_vtkActor SetVisibility $visible
    }

    method GetVisibility {} {
	return [$_vtkActor GetVisibility]
    }

    method GetProperty {} {
	return [$_vtkActor GetProperty]
    }

    method GetBounds {} {
	return [$_vtkActor GetBounds]
    }

    method GetActor {} {
	return $_vtkActor
    }

    method GetClassName {} {
	return "cvisActor"
    }

    method GetOrigin {} {
	return [$_vtkActor GetOrigin]
    }

    method SetOrigin {x y z} {
	$_vtkActor SetOrigin $x $y $z
    }

    method GetPosition {} {
	return [$_vtkActor GetPosition]
    }

    method SetPosition {x y z} {
	$_vtkActor SetPosition $x $y $z
    }

    method SetPickable {pickable} {
	$_vtkActor SetPickable $pickable
    }
}

class cvisActor {
    inherit cvisActorBase
    constructor {} {
    } {
	set _vtkActor [vtkActor $this.actor]
    }
}

class cvisLODActor {
    inherit cvisActorBase
    constructor {} {
    } {
	set _vtkActor [vtkLODActor $this.actor]
    }
}


class cvisWarpActorBase {
    inherit cvisActorBase

    variable ScaleFactor
    variable TotalScale

    constructor {} {
    } {
	foreach plane {x y z} {
	    set ScaleFactor($plane) 1.0
	    set TotalScale($plane) 1.0
	}
    }

    method ComputeScale {} {
	foreach plane {x y z} {
	    set TotalScale($plane) [expr $ScaleFactor($plane) * $Scale($plane)]
	}
	$_vtkActor SetScale $TotalScale(x) $TotalScale(y) $TotalScale(z)
    }

    method SetScaleFactor {x y z} {
	set ScaleFactor(x) $x
	set ScaleFactor(y) $y 
	set ScaleFactor(z) $z
	ComputeScale
    }

    method SetScale {x y z} {
	set Scale(x) $x
	set Scale(y) $y 
	set Scale(z) $z
	ComputeScale
    }
}

class cvisWarpActor {
    inherit cvisWarpActorBase
    constructor {} {
    } {
	set _vtkActor [vtkActor $this.actor]
    }
}

class cvisLODWarpActor {
    inherit cvisWarpActorBase
    constructor {} {
    } {
	set _vtkActor [vtkLODActor $this.actor]
    }
}


