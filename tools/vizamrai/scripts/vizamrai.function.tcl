##############################################################################
# Variable names to process, perhaps these should be read in from the 
# command line?
##############################################################################


##############################################################################
# Disable window from tk
##############################################################################
if {[winfo exists .]} {
    wm withdraw .
}

set variable_titles {Density}

##############################################################################
# Construct non-whitespace version of names for use as TCL objects
##############################################################################
set variable_names {}
foreach variable $variable_titles {
    regsub -all " " $variable "" result
    lappend variable_names $result
}

##############################################################################
# Load in packages required by Vizamrai
##############################################################################
package require vizamrai

##############################################################################
# Parse the command line arguments
# Should have better error checking here but hey this is a simple script
##############################################################################
# Get the filename from the command line
if {[file exists [lindex $argv 0]]} {
    set InputFileName [lindex $argv 0]
} {
    puts "Error: Can't read file <$parsed_args(filename)>"
    exit
}

set OutputFileName [lindex $argv 1]



##############################################################################
# The file reader objects, one for each variable
##############################################################################

puts "Computing the following variables:"
set index 0
foreach variable $variable_names {
    cvisReadSAMRAI $variable "$variable"
    $variable SetFileName $InputFileName
    puts "\t[lindex $variable_titles $index]"
    $variable SetCurrentVariableName [lindex $variable_titles $index]
    incr index
}

##############################################################################
# The slice extraction objects, one for each variable
##############################################################################
set SlicePlane "x"

foreach variable $variable_names {
    cvisExtractSlicePlane $variable.$SlicePlane "extractsliceplane.$SlicePlane"
    $variable.$SlicePlane SetSlicePlane x
    $variable.$SlicePlane SetInput $variable
}

# Read in the data.
# This is stuff for the Vizamrai work queue, Update makes objects
# do stuff and the work queue needs to be flushed.
foreach variable $variable_names {
    $variable Update
}
::cvis::queue::Execute

##############################################################################
# User Defined Function to compute something
##############################################################################
class myFunction {
    inherit cvisFunction
    constructor {name} {
	cvisFunction::constructor $name
    } {
    }

    variable VariableNames
    variable Results

    method SetVariableNames {names} {
	set VariableNames $names
    }

    method GetVariableNames {} {
	return $VariableNames
    }

    method SetResults {value} {
	foreach variable $VariableNames {
	    set Results($variable) $value
	}
    }

    method GetResults {variable} {
	return $Results($variable)
    }

    # Delta's are for the cell being intersected, args is an array 
    # of cell values for the variables the user wants to operate on.  See 
    # SetInput below on the function object below.

    # Arguments need to be in the same order as they are listed in 
    # VariableNames
    method Function {delta_x delta_y delta_z args} {

	set index 0
	foreach variable $VariableNames {
	    set value [lindex $args $index]
	    set Results($variable) [expr $Results($variable) + \
		    ($delta_y * $delta_z * $value)]
	    incr index
	}

	return
    }
}

# Instantiate an object of the type just created
myFunction function Function 

function SetVariableNames $variable_names

# Set the inputs, each cell value is passed via the args array with the 
# same index used here.  
set index 0
foreach variable $variable_names {
    function SetInput $index $variable.$SlicePlane
    incr index
}


##############################################################################
#  Some coordinates/spacing information
##############################################################################
# FYI There are three sets of coordinates floating around in 
# Vizamrai:
# screen - the display coordinate space for VTK/OpenGL, usually we work in 
#          this space for Vizamrai objects
# real   - the users problem space
# index  - this is the finest level grid spacing, used to help
#          deal with round off problems etc.
# 
# Can't use real as screen since some OpenGL implementations have
# some problems with large/small coordinate space numbers.

# Conversion from Real to Screen coordinates scaling factor
set Scaling [Density GetScaling]

# The bounds in finest grid index space.
scan [Density GetIndexBounds] "%d %d %d %d %d %d" \
	IndexBounds(min,x) IndexBounds(max,x) \
	IndexBounds(min,y) IndexBounds(max,y) \
	IndexBounds(min,z) IndexBounds(max,z)

# Spacing for the finest grid.
scan [Density GetIndexScaling] "%f %f %f" \
	Delta(x) \
	Delta(y) \
	Delta(z) 

# Get bounds in screen coordinates and convert to real space
scan [Density GetBounds] "%f %f %f %f %f %f" \
	RealBounds(min,x) RealBounds(max,x) \
	RealBounds(min,y) RealBounds(max,y) \
	RealBounds(min,z) RealBounds(max,z)

foreach plane {x y z} {
    set RealBounds(min,$plane) [expr $RealBounds(min,$plane) / $Scaling]
    set RealBounds(max,$plane) [expr $RealBounds(max,$plane) / $Scaling]
}


##############################################################################
# Loop over all the slices along X, writing the computed value to the file
##############################################################################

# Open the output file 
set outputfile [open $OutputFileName "w"]

# Loop over the finest grid index space
for {set index $IndexBounds(min,x)} {$index < $IndexBounds(max,x)} {incr index} {

    # Compute real and screen coordinates for this slice (middle of the cell)
    set real [expr (($index + 0.5) * $Delta(x)) \
	    + $RealBounds(min,x) ]
    set screen [expr $real * $Scaling]

    # Initialize the accumators in our function
    function SetResults 0.0

    # Update the slices (slice is set using in screen coordinates)
    foreach variable $variable_names {
	$variable.$SlicePlane SetSlice $screen
    }
    ::cvis::queue::Execute

    # Force computation of the users function
    function Update

    puts -nonewline $outputfile "$real "
    foreach variable $variable_names {
	set result [function GetResults $variable]
	set result [expr $Delta(x) * $result]

	set volume [expr ($RealBounds(max,y) - $RealBounds(min,y)) * \
		($RealBounds(max,z) - $RealBounds(min,z)) * \
		$Delta(x)]
	set result [expr $result / $volume]
    
	# Output the result
	puts -nonewline $outputfile "$result "
    }
    puts $outputfile ""
}

close $outputfile

exit



