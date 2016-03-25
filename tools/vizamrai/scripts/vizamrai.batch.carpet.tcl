#=============================================================================
# Load in the classes need by vizamrai
#=============================================================================
package require vizamrai 

# Color lookup table controls the mapping of values to colors
cvisLookupTable colorTable "ColorTable"
# This sets the color range from red to blue
colorTable SetHueRange 0.66667 0.0

#=============================================================================
# Parse input file which contains setting for this run
#=============================================================================
set filename [cvisParseArgs]

if {[file exists $filename]} {
    source $filename
} {
    puts "Error: Can't read file <$filename>"
    exit
}

#=============================================================================
# Create the VTK and cvis objects 
#=============================================================================

# Rendering object
cvisRenderWindow renderWindow "RenderWindow"
renderWindow SetOffScreenRendering 1
eval renderWindow SetBackground $Background

# This is colorbar which is drawn in the image to show
# what color to value mapping is
cvisColorBar colorBar "ColorBar"
colorBar SetLookupTable colorTable
colorBar SetRenderWindow renderWindow
colorBar SetLookupTable colorTable
eval colorBar SetColor $ColorBarTextColor

# Object which reads in SAMRAI data
cvisReadSAMRAI readSamrai "readSamrai"
readSamrai AddDepends colorBar
# Turn off autoscaling and use a constant value for min max data
# to cap data
eval readSamrai SetScalarRange $ScalarRange
readSamrai SetAutoScaleRange 0

# Surface plot object
cvisCarpetSlice carpetSlice "carpetSlice"
carpetSlice SetInput readSamrai
carpetSlice AddDepends renderWindow
carpetSlice SetLookupTable colorTable
carpetSlice SetRenderWindow renderWindow
carpetSlice SetBoundaryEdgeWidth $CarpetBoundaryEdgeWidth
carpetSlice SetCellToNode $CarpetCellToNode

cvisBoundingBox boundingBox "BoundingBox"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai
boundingBox SetBoundaryEdgeWidth $BoundingBoxBoundaryEdgeWidth

# Draw the axis indicator
#cvisAxes axes Axes
#axes SetRenderWindow renderWindow
#axes AddInput readSamrai

#=============================================================================
# Set the rendering view and image size
#=============================================================================

# Restore the camera position
renderWindow RestoreView $CameraFile
# Size of the image to create
eval renderWindow SetSize $ImageSize

#axes Update

#=============================================================================
# Loop over all the images in the supplied pattern
# and render each image to a file.
#=============================================================================
puts -nonewline "Working on frame "
set frame 0
# glob does file pattern matching
foreach filename [lsort [glob $FileNamePattern]] {
    puts -nonewline "..$frame"
    flush stdout
    if {[string length $filename]} {
	if {[file exists $filename]} {
	    #=================================================================
	    # Read the data from the file
	    #=================================================================
	    readSamrai SetFileName $filename
	    readSamrai Update

	    #=================================================================
	    # These two things need only be done on 1st dataset
	    # but can't be done until the dataset is loaded so we 
	    # know what the domain extent is and colormap needs 
	    # to be loaded by the dataset for the colorBar to work
	    # correctly
	    #=================================================================
	    if { $frame == 0 } {
		carpetSlice SetSlice $SlicePlane
		carpetSlice SetScaleFactor $CarpetScaleFactor
		colorBar Update
	    }
	    #=================================================================
	    # Render the image and save the image to a file
	    #=================================================================
	    renderWindow Render
	    renderWindow SetFileName "[file rootname $filename].ppm"
	    renderWindow SaveImageAsPPM
	} {
	    puts "Error: Can't read file <$filename>"
	    exit
	}

	incr frame
    }
}
puts ""
puts "done"

exit


