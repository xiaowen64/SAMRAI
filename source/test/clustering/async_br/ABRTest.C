/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2002 The Regents of the University of California
  Revision:	$Revision: 5 $
  Modified:	$Date: 2004-11-30 08:48:55 -0800 (Tue, 30 Nov 2004) $
  Description:	ABRTest class implementation
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ABRTest.h"
#include "CartesianGridGeometry.h"
#include "HierarchyCellDataOpsReal.h"
#include "ArrayData.h"
#include "CellVariable.h"
#include "NodeData.h"
#include "NodeVariable.h"

#include <iomanip>

using namespace SAMRAI;

template<int DIM>
ABRTest<DIM>::ABRTest(
  const string &object_name
, tbox::Pointer<tbox::Database> database
)
: d_name(object_name) ,
  d_hierarchy() ,
  d_tagger( object_name+":tagger",
	    database->isDatabase("sine_tagger") ?
	    database->getDatabase("sine_tagger").getPointer() : NULL ) ,
  d_time(0.5)
{
  return;
}

template<int DIM>
ABRTest<DIM>::~ABRTest()
{
  return;
}



template<int DIM>
mesh::StandardTagAndInitStrategy<DIM> *ABRTest<DIM>::getStandardTagAndInitObject()
{
  return &d_tagger;
}




/*
  Deallocate patch data allocated by this class.
*/
template<int DIM>
void ABRTest<DIM>::computeHierarchyData(
   hier::PatchHierarchy<DIM> &hierarchy,
   double time )
{
  d_tagger.computeHierarchyData( hierarchy, time );
  return;
}




/*
  Deallocate patch data allocated by this class.
*/
template<int DIM>
void ABRTest<DIM>::deallocatePatchData(
   hier::PatchHierarchy<DIM> &hierarchy )
{
  d_tagger.deallocatePatchData( hierarchy );
  return;
}





/*
  Deallocate patch data allocated by this class.
*/
template<int DIM>
void ABRTest<DIM>::deallocatePatchData(
   hier::PatchLevel<DIM> &level )
{
  d_tagger.deallocatePatchData( level );
  return;
}




template<int DIM>
int ABRTest<DIM>::registerVariablesWithPlotter(
  tbox::Pointer<appu::CartesianVizamraiDataWriter<DIM> > writer )
{
  if ( ! writer.isNull() )
    d_tagger.registerVariablesWithPlotter(*writer);
  return 0;
}




template<int DIM>
int ABRTest<DIM>::registerVariablesWithPlotter(
  tbox::Pointer<appu::VisItDataWriter<DIM> > writer )
{
  if ( ! writer.isNull() )
    d_tagger.registerVariablesWithPlotter(*writer);
  return 0;
}


template<int DIM>
bool ABRTest<DIM>::packDerivedDataIntoDoubleBuffer(
  double* buffer ,
  const hier::Patch<DIM> &patch ,
  const hier::Box<DIM> &region ,
  const string &variable_name ,
  int depth_id)
{
  if ( variable_name == "Patch<DIM> level number" ) {
    double pln = patch.getPatchLevelNumber();
    int i, size=region.size();
    for ( i=0; i<size; ++i ) buffer[i] = pln;
  }
  else {
    // Did not register this name.
    TBOX_ERROR("Unregistered variable name '" << variable_name << "' in\n"
	       << "ABRTest<DIM>::packDerivedPatchDataIntoDoubleBuffer");
  }

  return false;
}


#ifdef NDIM
template class ABRTest<NDIM>;
#endif
