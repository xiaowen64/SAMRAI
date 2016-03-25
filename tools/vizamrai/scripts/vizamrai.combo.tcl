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

cvisLookupTable colorTable "ColorTable"

cvisLookupTableInterface colorTableInterface "ColorTable"
colorTableInterface AddInput colorTable
colorTableInterface SetScalarRange 1e-6 1
colorTableInterface SetRenderWindow renderWindow


colorTable SetRenderWindow renderWindow

cvisRenderWindowInterface renderWindowInterface "RenderWindowInterface"
renderWindowInterface SetRenderWindow renderWindow

cvisReadSAMRAIInterface readSamraiInterface "ReadSAMRAIInterface"
readSamraiInterface SetRenderWindow renderWindow

mainWindow AddInterface readSamraiInterface
mainWindow AddInterface renderWindowInterface
mainWindow AddInterface colorTableInterface

cvisReadSAMRAI readSamrai "readSamrai"
readSamraiInterface SetInput readSamrai

cvisBoundingBoxInterface boundingBoxInterface "PatchBoundaries"
mainWindow AddInterface boundingBoxInterface 

cvisBoundingBox boundingBox "patchboundaries"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

#  cvisBoundingSlice boundingslice "BoundingSlice"
#  boundingslice SetRenderWindow renderWindow
#  boundingslice SetInput readSamrai

boundingBoxInterface AddInput boundingBox
boundingBoxInterface SetRenderWindow renderWindow

foreach plane "x y z" {
    cvisSlicePlaneInterface slicePlaneInterface.$plane \
	    "SlicePlaneInterface.$plane"
    mainWindow AddInterface slicePlaneInterface.$plane

    cvisSlicePlane slicePlane.$plane "SlicePlane.$plane"

    slicePlaneInterface.$plane AddInput slicePlane.$plane
    slicePlaneInterface.$plane SetSlicePlane $plane

    slicePlane.$plane SetInput readSamrai

    slicePlane.$plane AddDepends renderWindow
    
    slicePlane.$plane SetLookupTable colorTable
    
    slicePlane.$plane SetRenderWindow renderWindow
    slicePlaneInterface.$plane SetRenderWindow renderWindow
}

cvisVolInterface volInterface "Volume"
volInterface SetRenderWindow renderWindow
mainWindow AddInterface volInterface


cvisVol vol "vol"

vol SetInput readSamrai
vol SetRenderWindow renderWindow
vol SetLookupTable colorTable

volInterface AddInput vol
vol AddDepends renderWindow

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

colorTable SetInput readSamrai

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

# This is a hack to get the color table to redraw
# after some other geometry is on the screen
# otherwise it's color is not drawn correctly
colorTable SetVisibility 0
colorTable SetVisibility 1


::cvis::WaitForExit

foreach plane "x y z" {
    foreach i { slicePlane slicePlaneInterface} {
	::itcl::delete object $i.$plane
    }
}

foreach i { vol volInterface browser browserInterface axes axesInterface boundingBox boundingBoxInterface readSamrai readSamraiInterface renderWindowInterface colorTable colorTableInterface mainWindow renderWindow} {
    ::itcl::delete object $i
}

exit





