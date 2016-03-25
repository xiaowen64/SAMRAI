##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/Wrapping/Tcl/vizamrai/pkgIndex.tcl $
## Package:     Vizamrai
## Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Package loader for Vizamrai classes
##
##

# The cvis library path
# Is there a better place to put this?

global cvis_library
set cvis_library $dir

global cvisName

global cvisRelease
set cvisRelease [string trim "\$LastChangedRevision: 1704 $" "\$"]

global cvisDate
set cvisDate [string trim "\$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $" "\$"]

package ifneeded vizamrai 2.0 [list source [file join $dir vizamrai.tcl]]

# The cvis classes
package ifneeded cvis 2.0 [list tclPkgSetup $dir cvis 2.0 { \
	{cvisActor.tcl source {cvisActor cvisWarpActor}} \
	{cvisBaseClasses.tcl source {cvisDisplayObject cvisInterface}} \
	{cvisCollection.tcl source {cvisCollection}} \
	{cvisOperator.tcl source {cvisOperator}} \
	{cvisFunction.tcl source {cvisFunction}} \
	{cvisOperatorInterface.tcl source {cvisOperatorInterface}} \
	{cvisMain.tcl source {cvisMain}} \
	{cvisLookupTable.tcl source {cvisLookupTable}} \
	{cvisLookupTableInterface.tcl source {cvisLookupTableInterface}} \
	{cvisTabbedNotebook.tcl source {cvisTabbedNotebook}} \
	{cvisRenderWindowInterface.tcl source {cvisRenderWindowInterface}} \
	{cvisRenderWindow.tcl source {cvisRenderWindow}} \
	{cvis2DRenderWindow.tcl source {cvis2DRenderWindow}} \
	{cvisReadSAMRAIInterface.tcl source {cvisReadSAMRAIInterface}} \
	{cvisReadSAMRAI.tcl source {cvisReadSAMRAI}} \
	{cvisReadSAMRAIVectorInterface.tcl source {cvisReadSAMRAIVectorInterface}} \
	{cvisReadSAMRAIVector.tcl source {cvisReadSAMRAIVector}} \
	{cvisSlicePlaneInterface.tcl source {cvisSlicePlaneInterface}} \
	{cvisSlicePlane.tcl source {cvisSlicePlane}} \
	{cvisVector.tcl source {cvisVector}} \
	{cvisVectorInterface.tcl source {cvisVectorInterface}} \
	{cvisContourPlaneInterface.tcl source {cvisContourPlaneInterface}} \
	{cvisContourPlane.tcl source {cvisContourPlane}} \
	{cvisExtractSlicePlane.tcl source {cvisExtractSlicePlane}} \
	{cvisDataBrowserInterface.tcl source {cvisDataBrowserInterface}} \
	{cvisDataBrowser.tcl source {cvisDataBrowser}} \
	{cvisCarpetSliceInterface.tcl source {cvisCarpetSliceInterface}} \
	{cvisCarpetSlice.tcl source {cvisCarpetSlice}} \
	{cvisIsoInterface.tcl source {cvisIsoInterface}} \
	{cvisIso.tcl source {cvisIso}} \
	{cvisVol.tcl source {cvisVol}} \
	{cvisVolInterface.tcl source {cvisVolInterface}} \
	{cvisColorBar.tcl source {cvisColorBar}} \
	{cvisBoundingBox.tcl source {cvisBoundingBox}} \
	{cvisBoundingBoxInterface.tcl source {cvisBoundingBoxInterface}} \
	{cvisBoundingSlice.tcl source {cvisBoundingSlice}} \
	{cvisHackedBoundingBox.tcl source {cvisHackedBoundingBox}} \
	{cvisHackedBoundingBoxInterface.tcl source {cvisHackedBoundingBoxInterface}} \
	{cvisAxes.tcl source {cvisAxes}} \
	{cvisAxesInterface.tcl source {cvisAxesInterface}} \
	{cvisProperty.tcl source {cvisProperty}} \
	{cvisList.tcl source {cvisList}} \
	{cvisFileBrowser.tcl source {cvisFileBrowser}} \
	{cvisTextBrowser.tcl source {cvisTextBrowser}} \
	{cvisColorChooser.tcl source {cvisColorChooser}} \
	{cvisPlaneCellPatches.tcl source {cvisPlaneCellPatches}} \
	{foxTypedOpts.tcl source {typedopts}} \
	{cvisParseArgs.tcl source {cvisParseArgs}} \
	{cvisTrace.tcl source {cvisTrace}} \
	{cvisImageWriter.tcl source {cvisImageWriter}} \
    }]

# From the VTK package, put here so we can include it easier
package ifneeded vtkinteractor 1.0 [list tclPkgSetup $dir vtkinteractor 1.0 { \
    {vtkinteractor.tcl source {.vtkInteract .vtkInteract.buttons \
    .vtkInteract.buttons.dismiss .vtkInteract.display \
    .vtkInteract.display.scroll .vtkInteract.display.text .vtkInteract.file \
    .vtkInteract.file.entry .vtkInteract.file.label BindTkRenderWidget \
    EndMotion Enter Pan Render Reset Rotate StartMotion Surface \
    UpdateRenderer Wireframe Zoom dovtk vtkInteract}}}]






