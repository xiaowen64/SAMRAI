option add *textBackground white

set Bounds ""

# Load in the VTK classes
catch {load vtktcl}

# Load the ITCL package for OO capabilities
package require Itcl
namespace import itcl::*

# Advanced widgets from the [incr Widgets] package
package require Iwidgets 3.0

# tkcon package for nice command line interface
package require tkcon

# The cvis routines
package require cvis

# Parse the command line arguments
cvisParseArgs parsed_args

::cvis::queue::SetPipeline 0

cvisMain mainWindow ""
cvisRenderWindow renderWindow "RenderWindow"

cvisLookupTable colorTable "ColorTable"

cvisLookupTableInterface colorTableInterface "ColorTable"
colorTableInterface AddInput colorTable
colorTableInterface SetScalarRange 1e-6 1
colorTableInterface SetRenderWindow renderWindow

cvisRenderWindowInterface renderWindowInterface "RenderWindowInterface"
renderWindowInterface SetRenderWindow renderWindow

# Need to split interface and object
cvisReadSAMRAIInterface readSamraiInterface "ReadSAMRAIInterface"
readSamraiInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAI readSamrai "readSamrai"

readSamraiInterface SetInput readSamrai

colorTable SetInput readSamrai

cvisVolInterface volInterface "Volume"
volInterface SetRenderWindow renderWindow
mainWindow AddInterface volInterface

cvisVol vol "vol"

vol SetInput readSamrai
vol SetRenderWindow renderWindow
vol SetLookupTable colorTable

volInterface AddInput vol

vol AddDepends renderWindow

cvisBoundingBoxInterface boundingBoxInterface "BoundingBoxInterface"
mainWindow AddInterface boundingBoxInterface 

cvisBoundingBox boundingBox "BoundingBox"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

boundingBoxInterface AddInput boundingBox
boundingBoxInterface SetRenderWindow renderWindow

cvisAxes axes Axes
axes SetRenderWindow renderWindow
axes AddInput readSamrai
readSamrai AddDepends axes

cvisAxesInterface axesInterface AxesInterface
mainWindow AddInterface axesInterface 
axesInterface AddInput axes
axesInterface SetRenderWindow renderWindow

renderWindow SetSize 200 200

readSamrai SetEquivFileNameProp [renderWindow GetWindowNameProp]
readSamrai SetEquivFileNameProp [renderWindowInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [mainWindow GetWindowNameProp]

axes Update

mainWindow Pack

if {[string length [array names parsed_args filename]]} {
    if {[file exists $parsed_args(filename)]} {
	readSamraiInterface SetFileName $parsed_args(filename)
	set ScalarRange [readSamrai GetScalarRange]
	set mid [expr \
		( [lindex  $ScalarRange 0] + [lindex $ScalarRange 1] ) / 2.0]
    } {
	puts "Error: Can't read file <$parsed_args(filename)>"
	exit
    }
}

if {[string length [array names parsed_args restart]]} {
    # Load in the initial state from a save file
    if {[file exists $parsed_args(restart)]} {
	::cvis::RestartFile::Open $parsed_args(restart)
    } {
	puts "Error: Can't read file <$parsed_args(filename)>"
	exit
    }
} {
    # Force a window reset to get a reasonable initial view if
    # user is not restoring state from a file
    renderWindow BirdsEye
    renderWindow Render
}

# This is a hack to get the color table to redraw
# after some other geometry is on the screen
# otherwise it's color is not drawn correctly
colorTable SetVisibility 0
colorTable SetVisibility 1

::cvis::queue::SetPipeline 1
renderWindow Reset

colorTable SetTableRange 1 2
colorTable SetAutoTableRangeScaling 0
colorTable SetType 1

proc DoChange { } {
    set start 1
    set end   1e14
    set slices 10

    set spacing [expr log($end - $start) / $slices]

    for {set frame 0} {$frame < $slices} {incr frame 1} {
	set min_range [expr $start + exp($frame * $spacing)]
	set max_range [expr $start + exp(($frame+4) * $spacing)]

	colorTable SetTableRange $min_range $end
	renderWindow Render
    }
}

puts "Total Rendering Time = [time DoChange]"

exit