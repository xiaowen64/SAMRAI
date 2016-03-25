option add *textBackground white

# Load in packages required by vizamria
catch {load vtktcl}
package require Itcl
namespace import itcl::*
package require Iwidgets 3.0
package require tkcon
package require cvis

# Parse the command line arguments
cvisParseArgs parsed_args

# This can be turned on to trace the 
# execution; for debugging only
#cvisTrace::ACTIVATE all

option add *textBackground white

::cvis::queue::SetPipeline 0

cvisMain mainWindow ""

cvisRenderWindow renderWindow "RenderWindow"

cvisRenderWindowInterface renderWindowInterface "RenderWindowInterface"
renderWindowInterface SetRenderWindow renderWindow

cvisLookupTable colorTable "ColorTable"

cvisLookupTableInterface colorTableInterface "ColorTable"
colorTableInterface AddInput colorTable
colorTableInterface SetRenderWindow renderWindow

colorTable SetRenderWindow renderWindow

# Need to split interface and object
cvisReadSAMRAIInterface readSamraiInterface "ReadSAMRAIInterface"
readSamraiInterface SetRenderWindow renderWindow

cvisCarpetSliceInterface carpetSliceInterface "carpetSliceInterface"
carpetSliceInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface carpetSliceInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAI readSamrai "readSamrai"

readSamraiInterface SetInput readSamrai

cvisCarpetSlice carpetSlice "carpetSlice"


colorTable SetInput readSamrai

carpetSliceInterface AddInput carpetSlice

carpetSlice SetInput readSamrai
carpetSlice AddDepends renderWindow
carpetSlice SetLookupTable colorTable

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

carpetSlice SetRenderWindow renderWindow

cvisAxesInterface axesInterface AxesInterface
mainWindow AddInterface axesInterface 
axesInterface AddInput axes
axesInterface SetRenderWindow renderWindow

renderWindow SetSize 640 480

readSamrai SetEquivFileNameProp [renderWindow GetWindowNameProp]
readSamrai SetEquivFileNameProp [renderWindowInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [mainWindow GetWindowNameProp]

axes Update

mainWindow Pack



if {[string length [array names parsed_args filename]]} {
    if {[file exists $parsed_args(filename)]} {
	readSamraiInterface SetFileName $parsed_args(filename)
	carpetSliceInterface SetSlice 0
    } {
	puts "Error: Can't read file <$parsed_args(filename)>"
	exit
    }
}

# Restart the rendering pipeline
::cvis::queue::SetPipeline 1

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
    renderWindow Reset
    renderWindow Render
}

# This is a hack to get the color table to redraw
# after some other geometry is on the screen
# otherwise it's color is not drawn correctly
colorTable SetVisibility 0
colorTable SetVisibility 1

renderWindow BirdsEye

scan [readSamrai GetBounds] "%f %f %f %f %f %f" \
	real(min,x) real(max,x) \
	real(min,y) real(max,y) \
	real(min,z) real(max,z)

set plane [expr ($real(max,z) - $real(min,z)) / 2.0]

carpetSliceInterface SetSlice $plane


proc DoScaleFactor { } {
    foreach scale {0.5 0.75 1.0 1.25 1.5} {
	carpetSliceInterface SetScaleFactor $scale
	renderWindow Render
    }
}

puts "Total Scaling Time = [time DoScaleFactor]"

exit

