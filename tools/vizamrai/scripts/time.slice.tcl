lappend auto_path u:/users/smithsg/VTK/vizamrai/cvis

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

cvisRenderWindow renderWindow "RenderWindow"

cvisLookupTable colorTable "ColorTable"

colorTable SetRenderWindow renderWindow

cvisReadSAMRAI readSamrai "readSamrai"

colorTable SetInput readSamrai

cvisBoundingBox boundingBox "BoundingBox"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

foreach plane "z" {

    cvisSlicePlane slicePlane.$plane "SlicePlane.$plane"

    slicePlane.$plane SetInput readSamrai

    slicePlane.$plane AddDepends renderWindow
    
    slicePlane.$plane SetLookupTable colorTable
    
    slicePlane.$plane SetRenderWindow renderWindow
}

cvisAxes axes Axes
axes SetRenderWindow renderWindow
axes AddInput readSamrai
readSamrai AddDepends axes

renderWindow SetSize 640 480

axes Update

mainWindow Pack

if {[string length [array names parsed_args filename]]} {
    if {[file exists $parsed_args(filename)]} {
	readSamrai SetFileName $parsed_args(filename)
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

proc DoSliceChange { } {

    scan [readSamrai GetBounds] "%f %f %f %f %f %f" \
		real(min,x) real(max,x) \
		real(min,y) real(max,y) \
		real(min,z) real(max,z)

    

    set plane $real(min,z) 
    while { $plane < $real(max,z)} {
	slicePlane.z SetSlice $plane
	renderWindow Render
	update
	set plane [expr $plane + ($real(max,z) - $real(min,z)) / 10]
    }
}

wm withdraw .
update 

puts "Total Slicing Time = [time DoSliceChange]"
