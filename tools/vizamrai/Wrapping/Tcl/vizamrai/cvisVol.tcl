##
## File:        cvisVol.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: Volume render
##


class cvisVol {
    inherit cvisDisplayObject
    
    variable Reader ""

    variable PreviousCollection ""
    
    variable Interface ""

    variable VolumeExists 0
    variable VolumeIsVisible 1

    variable TableRange "0.0001 1.0"
    variable LookupTable ""

    variable Cropping 0
    variable CroppingRegionPlanes

    variable Type 0

    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	set CroppingRegionPlanes(x,min) 0
	set CroppingRegionPlanes(x,max) 1

	set CroppingRegionPlanes(y,min) 0
	set CroppingRegionPlanes(y,max) 1

	set CroppingRegionPlanes(z,min) 0
	set CroppingRegionPlanes(z,max) 1
    }

    method SetInput { reader } {
	set Reader $reader
	$reader AddOutput $this
    }

    method SetInterface { interface } {
	set Interface $interface
	$interface Update
    }

    method SetRenderWindow {window} {
	set RenderWindow $window
    }

    method SetTableRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]
	set TableRange "$min $max"

	::cvis::queue::Add "$this UpdateLookupTable"
    }

    method GetTableRange {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $TableRange 
    }

    method SetLookupTable {lut} {
	$lut AddDepends $this
	$lut AddLogDepends $this
	set LookupTable $lut
	UpdateLookupTable
    }

    
    method UpdateLookupTable {} {
	Update
    }

    variable ColorFileName
    method SetColorFileName { filename } {
	set ColorFileName $filename

	if {[string length $filename]} { 
	    set fp [open $filename r]
	    if { $fp == "" } {
		puts "Error opening colormap file <$filename>"
	    } else {

		if {!$ColorMapExists} {	
		    set ColorMapExists 1
		}

		for {set index 0} {$index < 256} {incr index} {
		    gets $fp line
		    scan $line "%i %i %i" r g b
			    [expr $r/256.0] \
			    [expr $g/256.0] \
			    [expr $b/256.0] 
		}
	    }
	}
    }

    method GetColorFileName {} {
	return $ColorFileName
    }

    variable AlphaFileName
    method SetAlphaFileName { filename } {
	set AlphaFileName $filename


	if {[string length $filename]} { 
	    set fp [open $filename r]
	    if { $fp == "" } {
		puts "Error opening alpha map file <$filename>"
	    } else {
		if {!$AlphaMapExists} {	
		    set AlphaMapExists 1
		}

		for {set index 0} {$index < 256} {incr index} {
		    gets $fp line
		    scan $line "%i %i %i" alpha a2 a3
			    [expr $alpha/256.0]
		}
	    }
	}
    }

    method GetAlphaFileName {} {
	return $AlphaFileName
    }

    method SetType {type} {
	set tracervar [::itcl::local cvisTrace #auto]
	if {$Type != $type} {
	    set Type $type
	    Update
	}
    }

    method GetType {} {
	return $Type
    }

    variable Depends ""

    method AddDepends { depend } {
	lappend Depends $depend
    }


    variable ColorMapExists 0
    variable AlphaMapExists 0
    method ExecuteVolumeCreate {} {
	if {![string length $RenderWindow] } { 
	    return
	}

	if {$VolumeExists} {
	    return
	}

	set VolumeExists 1

	if {!$AlphaMapExists} {	
	    set AlphaMapExists 1
	}

	if {!$ColorMapExists} {	
	    set ColorMapExists 1
	}
	
	
	
	set min [lindex $TableRange 0]
	set max [lindex $TableRange 1]
	
	set collection [$Reader GetCollection]
	$collection InitTraversal
	set item [$collection GetNextItem]
	while {[string length $item]} {

	    foreach plane {x y z} i {0 1 2}  {
		set origin($plane) [lindex [$item GetOrigin] $i]
		set dim($plane)    [lindex [$item GetDimensions] $i]
		set spacing($plane)    [lindex [$item GetSpacing] $i]
	    }

	    


	    if {$Type} {
	    } {
	    }

	    if {$Type} {
	    } {
	    }
	    
	    if {$Type} {
			[expr (255.0 / (log($max) - log($min)))]
	    } {
			[expr (255.0 / ($max - $min))]
	    }
	    
	    
	    

	    if {$Cropping} {
	    } {
	    }

		    $CroppingRegionPlanes(x,min) $CroppingRegionPlanes(x,max) \
		    $CroppingRegionPlanes(y,min) $CroppingRegionPlanes(y,max) \
		    $CroppingRegionPlanes(z,min) $CroppingRegionPlanes(z,max)


	    
	    
	    
	    
	    set item [$collection GetNextItem]
	}
	
	$RenderWindow SetScaling
    }

    method SetCroppingRegionPlanes {min_x max_x min_y max_y min_z max_z} {

	set CroppingRegionPlanes(x,min) $min_x
	set CroppingRegionPlanes(x,max) $max_x

	set CroppingRegionPlanes(y,min) $min_y
	set CroppingRegionPlanes(y,max) $max_y

	set CroppingRegionPlanes(z,min) $min_z
	set CroppingRegionPlanes(z,max) $max_z

	Update
    }

    method SetCropping {cropping} {
	if {$Cropping != $cropping} {
	    set Cropping $cropping

	    if {$VolumeExists} {
		set collection [$Reader GetCollection]
		$collection InitTraversal
		set item [$collection GetNextItem]
		
		while {[string length $item]} {
		    if {$Cropping} {
		    } {
		    }
		    set item [$collection GetNextItem]
		}
	    }
	}
    }

    method GetCropping {} {
	return $Cropping
    }

    method GetCroppingRegionPlanes {} {
	return "$CroppingRegionPlanes(x,min) \
		$CroppingRegionPlanes(x,max) \
		$CroppingRegionPlanes(y,min) \
		$CroppingRegionPlanes(y,max) \
		$CroppingRegionPlanes(z,min) \
		$CroppingRegionPlanes(z,max)" 
    }

    method ExecuteVolumeDelete {} {
	if {[string length \
	}
	    	    
	foreach item $PreviousCollection {
	    if {[string length \

	    }
	}
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]

	# Update user interface elements
	if {[string length $Interface]} {
	    $Interface Update
	}

	if {![string length $RenderWindow] } { 
	    return
	}

	# Schedule CPU intensive tasks that accomplish the update
	if {$VolumeExists} {
	    set VolumeExists 0
	    ::cvis::queue::Add "$this ExecuteVolumeDelete"
	}

	# Reset the collection information
	::cvis::queue::Add "$this ExecuteSetPreviousCollection"

	if {$VolumeIsVisible} {
	    ::cvis::queue::Add "$this ExecuteVolumeCreate"
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

}


