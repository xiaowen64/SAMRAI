##
## File:        cvisReadSAMRAI.tcl
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
## Description: File Reader for Vizamrai data
##

class cvisReadSAMRAI { 
    inherit cvisDisplayObject

    variable FileNameProp ""
    variable TimeProp ""

    variable VariableName ""
    variable VariableNames ""

    variable Reader ""

    variable Outputs ""

    variable Bounds "0 0 0 1 1 1"

    variable Depends ""

    variable ScalarRange "0 1"

    variable AutoScaleRange 1

    variable Cropping 0
    variable CroppingRegionPlanes



    method AddDepends { depend } {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Depends $depend
    }

    method GetBounds {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Bounds
    }

    method GetIndexBounds { } { 
	set tracervar [::itcl::local cvisTrace #auto]
	if {[string length $Reader]} {
	    return [$Reader GetIndexBounds]
	} {
	    return {0 1 0 1 0 1}
	}
    }

    method GetIndexScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetIndexScaling]
    }

    method SetStride {x y z} {
	$Reader SetStride $x $y $z
    }

    method GetStride {} {
	return [$Reader GetStride]
    }
    
    constructor {name} {
	cvisDisplayObject::constructor $name
    } {
	vtkSamraiStructuredPointsReader $name.cvisReadSAMRAI
	set Reader $name.cvisReadSAMRAI

	set FileNameProp [cvisProperty $name.FileNameProp]
	set TimeProp [cvisProperty $name.TimeProp]

	set CroppingRegionPlanes(x,min) 0
	set CroppingRegionPlanes(x,max) 1

	set CroppingRegionPlanes(y,min) 0
	set CroppingRegionPlanes(y,max) 1

	set CroppingRegionPlanes(z,min) 0
	set CroppingRegionPlanes(z,max) 1
   } 

   destructor {
   }

    method GetCollection { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetCollection]
    }

    method SetScalarRange {min max} {
	set tracervar [::itcl::local cvisTrace #auto]
	set ScalarRange "$min $max"
	Update
    }

    method SetAutoScaleRange { autoscalerange} {
	set tracervar [::itcl::local cvisTrace #auto]
	set AutoScaleRange $autoscalerange
	Update
    }

    method GetScalarRange { } {
	set tracervar [::itcl::local cvisTrace #auto]
	if { $AutoScaleRange } {
	    return [$Reader GetScalarRange]
	} {
	    return $ScalarRange
	}
    }

    method GetScaling {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetScaling]
    }

    method GetReader { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return $Reader
    }

    method AddOutput {list} {
	set tracervar [::itcl::local cvisTrace #auto]
	lappend Outputs $list
    }

    method GetTime {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$TimeProp GetValue]
    }

    method SetTime {time} {
	set tracervar [::itcl::local cvisTrace #auto]
	$TimeProp SetValue $time
    }

    method SetEquivTimeProp { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$TimeProp SetEquiv $property
    }

    method GetFileName { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$FileNameProp GetValue]
    }

    method SetEquivFileNameProp { property } {
	set tracervar [::itcl::local cvisTrace #auto]
	$FileNameProp SetEquiv $property
    }

    method GetNumberOfVariables { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetNumberOfVariables]
    }

    method GetNumberOfPatchBoundaries { } {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetNumberOfPatchBoundaries]
    }

    method GetPatchBoundary {i} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetPatchBoundary $i]
    }

    method GetPatchBoundaryLevel {i} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetPatchBoundaryLevel $i]
    }

    method GetVariableName {i} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetVariableName $i]
    }

    method GetDepth {i} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetDepth $i]
    }

    method SetIndex {index} {
	set tracervar [::itcl::local cvisTrace #auto]
	$Reader SetIndex $index
	Update
    }

    method GetIndex {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetIndex]
    }

    method SetCurrentVariableName {variablename} {
	set tracervar [::itcl::local cvisTrace #auto]
	$Reader SetCurrentVariableName $variablename
	Update
	return
    }

    method GetPatch {patch} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetPatch $patch]
    }

    method GetPatchLevel {i} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetPatchLevel $i]
    }

    method GetCurrentVariableName {} {
	set tracervar [::itcl::local cvisTrace #auto]
	return [$Reader GetCurrentVariableName]
    }

    method SetFileName { filename } {
	set tracervar [::itcl::local cvisTrace #auto]
	cvisTrace::print COMM "Reading $filename"

	$FileNameProp SetValue $filename

	$Reader SetFileName $filename

	# SGS what should we do if not filename

	if {[string length $filename]} { 
	    
	    # Read in the new data
	    $Reader Execute
	    
	    # Need to set the variable names from the file
	    # Need to set default variable name
	    set NumVariables [$Reader GetNumberOfVariables]
	    set i 0
	    for {set i 0} { $i < $NumVariables } {incr i} {
		lappend VariableNames [$Reader GetVariableName $i]
	    }
	    
	    # This does not seem to be always correctly updating the
	    # picked option.  Why is that?
	    set VariableName [lindex $VariableNames 0]

	    set Bounds [$Reader GetBounds]

	    SetTime [$Reader GetTime]
	}
    }

    method SetCroppingRegionPlanes {min_x max_x min_y max_y min_z max_z} {

	set CroppingRegionPlanes(x,min) $min_x
	set CroppingRegionPlanes(x,max) $max_x

	set CroppingRegionPlanes(y,min) $min_y
	set CroppingRegionPlanes(y,max) $max_y

	set CroppingRegionPlanes(z,min) $min_z
	set CroppingRegionPlanes(z,max) $max_z

	$Reader SetCroppingRegionPlanes \
		$CroppingRegionPlanes(x,min) $CroppingRegionPlanes(x,max) \
		$CroppingRegionPlanes(y,min) $CroppingRegionPlanes(y,max) \
		$CroppingRegionPlanes(z,min) $CroppingRegionPlanes(z,max) 
	
    }

    method SetCropping {cropping} {
	if {$Cropping != $cropping} {
	    set Cropping $cropping
	    $Reader SetCropping $cropping
	}
    }

    method Update {} {
	set tracervar [::itcl::local cvisTrace #auto]
	# Update things that depend on this
	foreach item $Depends {
	    $item Update
	}

	# Update children of this module
	foreach i $Outputs { 
	    $i SetInitialize
	    $i Update
	}
    }
}
