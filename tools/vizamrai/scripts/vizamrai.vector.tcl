package require vizamrai

# Parse the command line arguments
cvisParseArgs parsed_args

# This can be turned on to trace the 
# execution; for debugging only
# cvisTrace::ACTIVATE all

option add *textBackground white

::cvis::queue::SetPipeline 0

cvisMain mainWindow ""
cvisRenderWindow renderWindow "RenderWindow"

cvisLookupTable colorTable "colortable"

cvisLookupTableInterface colorTableInterface "ColorTable"
colorTableInterface AddInput colorTable
colorTableInterface SetRenderWindow renderWindow

colorTable SetRenderWindow renderWindow

cvisRenderWindowInterface renderWindowInterface "RenderWindow"
renderWindowInterface SetRenderWindow renderWindow

cvisReadSAMRAIVectorInterface readSamraiInterface "ReadSAMRAIVector"
readSamraiInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAIVector readSamrai "readsamrai"
readSamraiInterface SetInput readSamrai



cvisBoundingBoxInterface boundingBoxInterface "PatchBoundaries"
mainWindow AddInterface boundingBoxInterface 

cvisBoundingBox boundingBox "patchboundaries"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

boundingBoxInterface AddInput boundingBox
boundingBoxInterface SetRenderWindow renderWindow

cvisVectorInterface vectorInterface \
	"Vector"
mainWindow AddInterface vectorInterface

cvisVector vector "vector"

vectorInterface AddInput vector

vector SetInput readSamrai

vector AddDepends renderWindow

vector SetLookupTable colorTable

vector SetRenderWindow renderWindow
vectorInterface SetRenderWindow renderWindow

cvisAxes axes Axes
axes SetRenderWindow renderWindow
axes AddInput readSamrai
readSamrai AddDepends axes

cvisAxesInterface axesInterface AxesInterface
mainWindow AddInterface axesInterface 
axesInterface AddInput axes
axesInterface SetRenderWindow renderWindow


#colorTable SetInput readSamrai
colorTable SetInput vector

renderWindow SetSize 640 480

# Propogate filename changes
readSamrai SetEquivFileNameProp [renderWindow GetWindowNameProp]
readSamrai SetEquivFileNameProp [renderWindowInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [mainWindow GetWindowNameProp]

axes Update

mainWindow Pack

if {[string length [array names parsed_args filename]]} {
    if {[file exists $parsed_args(filename)]} {
	readSamraiInterface SetFileName $parsed_args(filename)
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

foreach i {vector vectorInterface} {
    ::itcl::delete object $i
}

foreach i { axes axesInterface boundingBox boundingBoxInterface readSamrai readSamraiInterface renderWindowInterface colorTable colorTableInterface mainWindow renderWindow} {
    ::itcl::delete object $i
}

exit





