option add *textBackground white

set Bounds ""

# Load in the VTK classes
package require vizamrai

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

cvisRenderWindowInterface renderWindowInterface "RenderWindow"
renderWindowInterface SetRenderWindow renderWindow

# Need to split interface and object
cvisReadSAMRAIInterface readSamraiInterface "ReadSAMRAI"
readSamraiInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAI readSamrai "readsamrai"

readSamraiInterface SetInput readSamrai

cvisVolInterface volInterface "Volume"
volInterface SetRenderWindow renderWindow
mainWindow AddInterface volInterface

cvisVol vol "vol"

vol SetInput readSamrai
vol SetRenderWindow renderWindow
vol SetLookupTable colorTable

volInterface AddInput vol

vol AddDepends renderWindow

cvisBoundingBoxInterface boundingBoxInterface "PatchBoundaries"
mainWindow AddInterface boundingBoxInterface 

cvisBoundingBox boundingBox "patchboundaries"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

boundingBoxInterface AddInput boundingBox
boundingBoxInterface SetRenderWindow renderWindow

cvisAxes axes "axes"
axes SetRenderWindow renderWindow
axes AddInput readSamrai
readSamrai AddDepends axes

cvisAxesInterface axesInterface "Axes"
mainWindow AddInterface axesInterface 
axesInterface AddInput axes
axesInterface SetRenderWindow renderWindow

renderWindow SetSize 200 200

readSamrai SetEquivFileNameProp [renderWindow GetWindowNameProp]
readSamrai SetEquivFileNameProp [renderWindowInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [mainWindow GetWindowNameProp]

axes Update

colorTable SetInput readSamrai

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

::cvis::WaitForExit

foreach i { vol volInterface axes axesInterface boundingBox boundingBoxInterface readSamrai readSamraiInterface renderWindowInterface colorTable colorTableInterface mainWindow renderWindow} {
    ::itcl::delete object $i
}

exit


