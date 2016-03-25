##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/cvisAxes.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Labeling for Axes
##


class cvisAxes {
    inherit cvisDisplayObject

    variable Input ""

    variable Bounds 
    
    variable AxesType 1

    variable Visible 1

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	foreach plane {x y z} {
	    set Bounds(min,$plane) 0
	    set Bounds(max,$plane) 1
	}
    }

    method SetAxesType {type} {
	if {$AxesType != $type} {
	    set AxesType $type
	}
	set Create 1
	Update
    }
    
    method SetVisibility {visible} {
	if {$Visible != $visible} {
	    set Visible $visible
	    if {$AxesType} {
		$this.actor SetVisibility $Visible
	    } {
	    }
	}
    }

    method GetAxesActor { } {
	set tracervar [::itcl::local cvisTrace #auto]

	if {[string length [info commands $this.actor]]} {
	    return $this.actor
	} {
	    return ""
	}
    }

    method AddInput {data} {
	set tracervar [::itcl::local cvisTrace #auto]


	lappend Input $data
	$data AddDepends $this

	# if no input then set the values to be what this input 
	# is
	Update
    }
    
    method Update {} { 
	set tracervar [::itcl::local cvisTrace #auto]

	# Compute the bounds for scaling values
	foreach plane "x y z" { 
	    # SGS where can we get these from, VTK class
	    set Bounds(min,$plane) $::cvis::math::FLT_MAX
	    set Bounds(max,$plane) [expr -$::cvis::math::FLT_MAX]
	}

	foreach item $Input {

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
	}

	# If the vtk classes are not created then make them
	if $Create { 
	    if {[string length [info commands $this.actor]]} {
		$RenderWindow RemoveActor $this.actor
		::itcl::delete object $this.actor 
	    }

	    }

	    if {$AxesType} {

		
		
		
		
		# SGS was LOD
		cvisActor $this.actor
		$this.actor SetVisibility $Visible
		$this.actor SetPickable 0
		
		$RenderWindow AddActor $this.actor
	    } {
		
			GetActiveCamera]
		
	    }
	}
	set Create 0

	# Set the scale based on what the input sizes where
	set scale 0
	foreach plane "x y z" {
	    set temp_scale [expr ($Bounds(max,$plane) - \
		    $Bounds(min,$plane)) * 0.1]
	    if {$scale < $temp_scale} {
		set scale $temp_scale
	    }
	}

	if {$AxesType} {

	} {
		    $Bounds(min,y) $Bounds(max,y) \
		    $Bounds(min,z) $Bounds(max,z) 
	}
	
	$RenderWindow SetScaling
    }
}
