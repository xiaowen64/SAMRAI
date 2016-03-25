##
## File:        cvisIso.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Isosurface 
##

class cvisIso {
    inherit cvisDisplayObject
    
    variable Reader ""

    variable IsoValue 0.5

    variable PreviousCollection ""
    
    variable TableRange "0.0 1.0"

    variable LookupTable ""

    variable Interface ""

    variable IsoSurfaceExists 0
    # SGS should add the ability to toggle this from interface
    variable IsoSurfaceIsVisible 1

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
    }

    method SetInput { reader } {
	set Reader $reader
	$reader AddOutput $this
    }

    method SetInterface { interface } {
	set Interface $interface
	$interface Update
    }

    method GetIsoValue {} {
	return $IsoValue
    }

    method ExecuteSetIsoValue {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {![string length $Reader]} { 
	    return
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	set collection [[$Reader GetReader] GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    
		
		if {![string length \
		    
		    
		    
			    $::cvis::options::LOWMEM
		    
		    
		    if {[string length $LookupTable] } { 
				[$LookupTable GetLookupTable]
		    }
		    
		    
		}
	    } {
		if {[string length \
		    
		}
	    }
	    
	    set item [$collection GetNextItem]
	}
    }

    method SetIsoValue {value} {
	if {$IsoValue == $value} {
	    return
	}

	set IsoValue $value
	::cvis::queue::Add "$this ExecuteSetIsoValue"
    }

    method SetRenderWindow {window} {
	set RenderWindow $window
    }

    method SetTableRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]

	set TableRange "$min $max"

	if {[string length $Interface]} {
	    $Interface SetTableRange $min $max
	}

	::cvis::queue::Add "$this ExecuteSetTableRange"
    }


    method ExecuteSetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]

    	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    }
	    set item [$collection GetNextItem]
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
	lappend Depends $depend
    }

    method ExecuteIsoSurfaceCreate {} {
	if {$IsoSurfaceExists} {
	    return
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	set ScalarRange [$Reader GetScalarRange]
	
	if { $IsoValue < [lindex $ScalarRange 0] } {
	    set IsoValue [lindex $ScalarRange 0]
	}

	if { $IsoValue > [lindex $ScalarRange 1] } {
	    set IsoValue [lindex $ScalarRange 1]
	}

	set IsoSurfaceExists 1
	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {
	    
	    
	    
	    
		
		
			$::cvis::options::LOWMEM
		
		if {[string length $LookupTable] } { 
			    [$LookupTable GetLookupTable]
		}
		
		
	    }
	    
	    set item [$collection GetNextItem]
	}
	
	$RenderWindow SetScaling

	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}
    }
    method ExecuteIsoSurfaceDelete {} {
	set IsoSurfaceExists 0	
	
	foreach item $PreviousCollection {
	    if {[string length \
	    }
	    if {[string length \

	    }
	}
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]

	if {![string length $RenderWindow] } { 
	    return
	}

	# Schedule CPU intensive tasks that accomplish the update
	if {$IsoSurfaceExists} {
	    set IsoSurfaceExists 0
	    ::cvis::queue::Add "$this ExecuteIsoSurfaceDelete"
	}

	# Reset the collection information
	::cvis::queue::Add "$this ExecuteSetPreviousCollection"

	if {$IsoSurfaceIsVisible} {
	    ::cvis::queue::Add "$this ExecuteIsoSurfaceCreate"
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




