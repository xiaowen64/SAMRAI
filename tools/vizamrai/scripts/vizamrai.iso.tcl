option add *textBackground white

# Load in packages required by vizamria
package require vizamrai

# Parse the command line arguments
cvisParseArgs parsed_args

# This can be turned on to trace the 
# execution; for debugging only
#cvisTrace::ACTIVATE all

::cvis::queue::SetPipeline 0

cvisMain mainWindow ""
cvisRenderWindow renderWindow "RenderWindow"

cvisLookupTable colorTable "ColorTable"
colorTable SetRenderWindow renderWindow

cvisLookupTableInterface colorTableInterface "ColorTable"
colorTableInterface AddInput colorTable
colorTableInterface SetRenderWindow renderWindow

cvisRenderWindowInterface renderWindowInterface "RenderWindow"
renderWindowInterface SetRenderWindow renderWindow

cvisReadSAMRAIInterface readSamraiInterface "ReadSAMRAI"
readSamraiInterface SetRenderWindow renderWindow

cvisIsoInterface isoInterface "Iso"
isoInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface isoInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAI readSamrai "readSamrai"

readSamraiInterface SetInput readSamrai

cvisIso iso "iso"

isoInterface AddInput iso

iso SetInput readSamrai
iso AddDepends renderWindow
iso SetLookupTable colorTable

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

iso SetRenderWindow renderWindow

renderWindow SetSize 640 480

readSamrai SetEquivFileNameProp [renderWindow GetWindowNameProp]
readSamrai SetEquivFileNameProp [renderWindowInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [mainWindow GetWindowNameProp]

colorTable SetInput readSamrai

axes Update

mainWindow Pack

if {[string length [array names parsed_args filename]]} {
    if {[file exists $parsed_args(filename)]} {
	readSamraiInterface SetFileName $parsed_args(filename)
	set ScalarRange [readSamrai GetScalarRange]
	set mid [expr \
		( [lindex  $ScalarRange 0] + [lindex $ScalarRange 1] ) / 2.0]
	iso SetIsoValue $mid
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

foreach i { iso isoInterface axes axesInterface boundingBox boundingBoxInterface readSamrai readSamraiInterface renderWindowInterface colorTable colorTableInterface mainWindow renderWindow} {
    ::itcl::delete object $i
}

exit



