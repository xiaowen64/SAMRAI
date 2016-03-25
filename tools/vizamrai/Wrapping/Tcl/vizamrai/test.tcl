#!/bin/sh
# the next line restarts using wish \
exec CVIS_INSTALL/bin/cvish "$0" "$@"

puts $auto_path

package require vtkcvis

vtkSamraiStructuredPointsReader foo

exit

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

colorTable SetInput readSamrai

cvisBoundingBoxInterface boundingBoxInterface "BoundingBox"
mainWindow AddInterface boundingBoxInterface 

cvisBoundingBox boundingBox "boundingbox"
boundingBox SetRenderWindow renderWindow
boundingBox SetInput readSamrai

boundingBoxInterface AddInput boundingBox
boundingBoxInterface SetRenderWindow renderWindow

foreach plane "z" {
    cvisExtractSlicePlane slicePlane.$plane "sliceplane.$plane"
    slicePlane.$plane SetSlicePlane $plane
    slicePlane.$plane SetInput readSamrai

}


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

slicePlane.$plane ExecuteDataCreate
cvisOperator operator Operator
operator SetFileName "aaa.foo"
operator SetInput 0 [slicePlane.$plane GetOutput]
operator Update

exit


