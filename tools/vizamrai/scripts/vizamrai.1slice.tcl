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

foreach plane "z" {
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

    cvisContourPlaneInterface contourPlaneInterface.$plane \
	    "ContourPlane.$plane"
    mainWindow AddInterface contourPlaneInterface.$plane

    cvisContourPlane contourPlane.$plane "contourplane.$plane"

    contourPlaneInterface.$plane AddInput contourPlane.$plane

    contourPlane.$plane SetInput readSamrai

    contourPlane.$plane AddDepends renderWindow
    
    contourPlane.$plane SetLookupTable colorTable
    
    contourPlane.$plane SetRenderWindow renderWindow
    contourPlaneInterface.$plane SetRenderWindow renderWindow

    slicePlaneInterface.$plane AddPlaneInput contourPlane.$plane
}

##############################################################################
# Extract Slice Section
##############################################################################
foreach plane "z" {

    cvisExtractSlicePlane extractSlicePlane.$plane "extractsliceplane.$plane"
    extractSlicePlane.$plane SetSlicePlane $plane
    extractSlicePlane.$plane SetInput readSamrai

    # Use interface from the displayed slicer to control the slice
    # plane and slice
    slicePlaneInterface.$plane AddInput extractSlicePlane.$plane

}

##############################################################################
# Operator section
##############################################################################

cvisOperatorInterface operatorInterface "Save Slice"
mainWindow AddInterface operatorInterface

class myOperator {
    inherit cvisOperator
    constructor {name} {
	cvisOperator::constructor $name
    } {
    }

    method Function {delta_x delta_y delta_z args} {
	set value [lindex $args 0]
#	set value [expr $value / ($delta_x * $delta_y * $delta_z)]
	return $value
    }
}

 

myOperator operator Operator
operator SetInput 0 extractSlicePlane.$plane

# When extract slice gets updated we should update operator
extractSlicePlane.$plane AddDepends operator

operatorInterface AddInput operator

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

colorTable SetInput readSamrai

renderWindow SetSize 640 480

# Propogate filename changes
readSamrai SetEquivFileNameProp [renderWindow GetWindowNameProp]
readSamrai SetEquivFileNameProp [renderWindowInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [operatorInterface GetFileNameProp]
readSamrai SetEquivFileNameProp [mainWindow GetWindowNameProp]

# Propogate time changes
readSamrai SetEquivTimeProp     [operator GetTimeProp]

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

foreach plane "z" {
    foreach i {contourPlane contourPlaneInterface slicePlane slicePlaneInterface} {
	::itcl::delete object $i.$plane
    }
}

foreach i { operator operatorInterface browser browserInterface axes axesInterface boundingBox boundingBoxInterface readSamrai readSamraiInterface renderWindowInterface colorTable colorTableInterface mainWindow renderWindow} {
    ::itcl::delete object $i
}

exit





