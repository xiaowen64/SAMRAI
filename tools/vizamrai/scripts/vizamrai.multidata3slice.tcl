# Load in packages required by vizamria
package require vizamrai

# Parse the command line arguments
cvisParseArgs parsed_args

# This can be turned on to trace the 
# execution; for debugging only
#cvisTrace::ACTIVATE all

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

cvisReadSAMRAIInterface readSamraiInterface "ReadSAMRAI"
readSamraiInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAI readSamrai "readsamrai"
readSamraiInterface SetInput readSamrai


cvisBoundingBoxInterface boundingBoxInterface "PatchBoundaries"
mainWindow AddInterface boundingBoxInterface 

cvisBoundingBox boundingBox "patchboundaries"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

boundingBoxInterface AddInput boundingBox
boundingBoxInterface SetRenderWindow renderWindow

foreach plane "x y" {
    cvisSlicePlaneInterface slicePlaneInterface.$plane \
	    "SlicePlane.$plane"
    mainWindow AddInterface slicePlaneInterface.$plane

    cvisSlicePlane slicePlane.$plane "sliceplane.$plane"

    slicePlaneInterface.$plane AddInput slicePlane.$plane
    slicePlaneInterface.$plane SetSlicePlane $plane

    slicePlane.$plane SetInput readSamrai

    slicePlane.$plane AddDepends renderWindow
    
    slicePlane.$plane SetLookupTable colorTable
    
    slicePlane.$plane SetRenderWindow renderWindow
    slicePlaneInterface.$plane SetRenderWindow renderWindow
}

cvisAxes axes Axes
axes SetRenderWindow renderWindow
axes AddInput readSamrai
readSamrai AddDepends axes

cvisAxesInterface axesInterface AxesInterface
mainWindow AddInterface axesInterface 
axesInterface AddInput axes
axesInterface SetRenderWindow renderWindow


cvisDataBrowser browser DataBrowser
browser SetRenderWindow renderWindow
browser SetInput readSamrai
browser AddDepends renderWindow

cvisDataBrowserInterface browserInterface DataBrowser
mainWindow AddInterface browserInterface
browserInterface AddInput browser
browserInterface SetRenderWindow renderWindow
renderWindow AddPicker browserInterface

renderWindow SetSize 640 480

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

# *************************************************************************
# Second dataset
# *************************************************************************
cvisLookupTable colorTable2 "colortable2"

cvisLookupTableInterface colorTableInterface2 "ColorTable2"
colorTableInterface2 AddInput colorTable2
colorTableInterface2 SetRenderWindow renderWindow


colorTable2 SetRenderWindow renderWindow

cvisReadSAMRAIInterface readSamraiInterface2 "ReadSAMRAI2"
readSamraiInterface2 SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface2
mainWindow AddInterface colorTableInterface2

cvisReadSAMRAI readSamrai2 "readsamrai2"
readSamraiInterface2 SetInput readSamrai2

foreach plane "z" {
    cvisSlicePlaneInterface slicePlaneInterface2.$plane \
	    "SlicePlane2.$plane"
    mainWindow AddInterface slicePlaneInterface2.$plane

    cvisSlicePlane slicePlane2.$plane "sliceplane2.$plane"

    slicePlaneInterface2.$plane AddInput slicePlane2.$plane
    slicePlaneInterface2.$plane SetSlicePlane $plane

    slicePlane2.$plane SetInput readSamrai2

    slicePlane2.$plane AddDepends renderWindow
    
    slicePlane2.$plane SetLookupTable colorTable2
    
    slicePlane2.$plane SetRenderWindow renderWindow
    slicePlaneInterface2.$plane SetRenderWindow renderWindow
}

colorTable2 SetInput readSamrai2
colorTable SetInput readSamrai


colorTableInterface2 SetVisibility 0

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







