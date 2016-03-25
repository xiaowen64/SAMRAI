/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dlbg/SinusoidalFrontTagger.C $
  Copyright:	(c) 1997-2002 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	SinusoidalFrontTagger class implementation
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "SinusoidalFrontTagger.h"
#include "MDA_Access.h"
#include "CartesianGridGeometry.h"
#include "HierarchyCellDataOpsReal.h"
#include "ArrayData.h"
#include "CellVariable.h"
#include "NodeData.h"
#include "NodeVariable.h"

#include <iomanip>

using namespace SAMRAI;



template<int DIM>
SinusoidalFrontTagger<DIM>::SinusoidalFrontTagger(
  const string &object_name
, tbox::Database *database
)
: d_name(object_name) ,
  d_hierarchy() ,
  d_period( 1.0 ) ,
  d_amplitude( 0.2 ) ,
  d_adaption_buffer(1),
  d_allocate_data(true) ,
  d_time(0.5)
{
  hier::VariableDatabase<DIM> *variable_db = hier::VariableDatabase<DIM>::getDatabase();
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT( variable_db != NULL );
#endif

  tbox::Array<double> init_disp;
  tbox::Array<double> velocity;

  if ( database != NULL ) {
    d_allocate_data
      = database->getBoolWithDefault( "allocate_data" ,
				      d_allocate_data );
    d_adaption_buffer
      = database->getIntegerWithDefault( "adaption_buffer" ,
					 d_adaption_buffer );
    d_period
      = database->getDoubleWithDefault( "period" ,
					 d_period );
    if ( database->isDouble("init_disp") ) {
       init_disp
          = database->getDoubleArray( "init_disp" );
    }
    if ( database->isDouble("velocity") ) {
       velocity
          = database->getDoubleArray( "velocity" );
    }
    d_amplitude
      = database->getDoubleWithDefault( "amplitude" ,
					 d_amplitude );
    d_time
      = database->getDoubleWithDefault( "time" ,
					 d_time );
  }

  for ( int idim=0; idim<DIM; ++idim ) {
     d_init_disp[idim] = idim < init_disp.size() ? init_disp[idim] : 0.0;
     d_velocity[idim] = idim < velocity.size() ? velocity[idim] : 0.0;
  }

  const string context_name = d_name + string(":context");
  d_context = variable_db->getContext(context_name);

  const hier::IntVector<DIM> ghost(0);

  tbox::Pointer<hier::Variable<DIM> > dist_var =
    new pdat::NodeVariable<DIM,double>(d_name+":dist");
  d_dist_id = variable_db->registerVariableAndContext( dist_var,
						       d_context,
						       ghost);

  tbox::Pointer<hier::Variable<DIM> > tag_var =
    new pdat::CellVariable<DIM,int>(d_name+":tag");
  d_tag_id = variable_db->registerVariableAndContext( tag_var,
						      d_context,
						      ghost);

  return;
}



template<int DIM>
SinusoidalFrontTagger<DIM>::~SinusoidalFrontTagger()
{
  return;
}



template<int DIM>
void SinusoidalFrontTagger<DIM>::initializeLevelData (
  /*! Hierarchy to initialize */
  const tbox::Pointer< hier::BasePatchHierarchy<DIM> > base_hierarchy ,
  /*! Level to initialize */
  const int ln ,
  const double init_data_time ,
  const bool can_be_refined ,
  /*! Whether level is being introduced for the first time */
  const bool initial_time ,
  /*! Level to copy data from */
  const tbox::Pointer< hier::BasePatchLevel<DIM> > old_base_level ,
  const bool allocate_data )
{

   tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy = base_hierarchy;
   tbox::Pointer<hier::PatchLevel<DIM> > old_level = old_base_level;
   if ( ! old_base_level.isNull() ) {
      TBOX_ASSERT( ! old_level.isNull() );
   }
   TBOX_ASSERT( ! hierarchy.isNull() );

  /*
    Reference the level object with the given index from the hierarchy.
  */
  tbox::Pointer<hier::PatchLevel<DIM> > level
    = hierarchy->getPatchLevel(ln);
  typename hier::PatchLevel<DIM>::Iterator pi;
  for ( pi.initialize(level); pi; pi++ ) {
     hier::Patch<DIM> &patch = *level->getPatch(*pi);
     initializePatchData( patch,
                          init_data_time,
                          initial_time,
                          allocate_data );
  }

#if 0
  if ( d_allocate_data ) {
    /*
      If instructed, allocate all patch data on the level.
      Allocate only persistent data.  Scratch data will
      generally be allocated and deallocated as needed.
    */
    if ( allocate_data ) {
      level->allocatePatchData( d_dist_id );
      level->allocatePatchData( d_tag_id );
    }
    computeLevelData( hierarchy, ln, d_time/*init_data_time*/,
		      d_dist_id, d_tag_id, old_level );
  }
#endif

  return;
}



template<int DIM>
void SinusoidalFrontTagger<DIM>::initializePatchData (
  /*! Hierarchy to initialize */
  hier::Patch<DIM> &patch ,
  const double init_data_time ,
  const bool initial_time ,
  const bool allocate_data )
{

  if ( d_allocate_data ) {
    /*
      If instructed, allocate all patch data on the level.
      Allocate only persistent data.  Scratch data will
      generally be allocated and deallocated as needed.
    */
    if ( allocate_data ) {
      patch.allocatePatchData( d_dist_id );
      patch.allocatePatchData( d_tag_id );
      tbox::Pointer<pdat::NodeData<DIM,double> > dist_data =
         patch.getPatchData(d_dist_id);
      tbox::Pointer<pdat::CellData<DIM,int> > tag_data =
         patch.getPatchData(d_tag_id);
      TBOX_ASSERT( ! dist_data.isNull() );
      TBOX_ASSERT( ! tag_data.isNull() );
      computePatchData( patch, init_data_time,
                        dist_data.getPointer(), tag_data.getPointer() );
    }
  }

  return;
}



template<int DIM>
void SinusoidalFrontTagger<DIM>::resetHierarchyConfiguration (
  /*! New hierarchy */ tbox::Pointer<hier::BasePatchHierarchy<DIM> > new_hierarchy ,
  /*! Coarsest level */ int coarsest_level ,
  /*! Finest level */ int finest_level )
{
  d_hierarchy = new_hierarchy;
  TBOX_ASSERT( ! d_hierarchy.isNull() );
  return;
}



template<int DIM>
void SinusoidalFrontTagger<DIM>::applyGradientDetector(
  const tbox::Pointer<hier::BasePatchHierarchy<DIM> > base_hierarchy_,
  const int ln,
  const double error_data_time,
  const int tag_index,
  const bool initial_time,
  const bool uses_richardson_extrapolation )
{
   tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy_ = base_hierarchy_;
   TBOX_ASSERT( ! hierarchy_.isNull() );
  tbox::Pointer<hier::PatchLevel<DIM> > level_ = hierarchy_->getPatchLevel(ln);
   TBOX_ASSERT( ! level_.isNull() );

  hier::PatchLevel<DIM> &level = *level_;

  typename hier::PatchLevel<DIM>::Iterator pi;
  for ( pi.initialize(level); pi; pi++ ) {
    const int pn = pi();
    hier::Patch<DIM> &patch = *level.getPatch(pn);

    tbox::Pointer<hier::PatchData<DIM>  >
      tag_data = patch.getPatchData( tag_index );
    if ( tag_data.isNull() ) {
      TBOX_ERROR("Data index " << tag_index
		 << " does not exist for patch.\n");
    }
    tbox::Pointer<pdat::CellData<DIM,int> > tag_cell_data_ = tag_data;
    if ( tag_cell_data_.isNull() ) {
      TBOX_ERROR("Data index " << tag_index
		 << " is not cell int data.\n");
    }

    if ( d_allocate_data ) {
      // Use internally stored data.
      tbox::Pointer<pdat::CellData<DIM,int> >
	saved_tag_data = patch.getPatchData( d_tag_id );
      tag_cell_data_->copy( *saved_tag_data );
      // overlayTagData( *saved_tag_data, *tag_cell_data_ );
    }
    else {
      // Compute tag data for patch.
      computePatchData( patch,
                        error_data_time,
                        NULL,
                        tag_cell_data_.getPointer() );
    }
	
  }

  return;
}




/*
  Deallocate patch data allocated by this class.
*/


template<int DIM>
void SinusoidalFrontTagger<DIM>::deallocatePatchData(
   hier::PatchHierarchy<DIM> &hierarchy )
{
  int ln;
  for ( ln=0; ln<hierarchy.getNumberLevels(); ++ln ) {
     tbox::Pointer<hier::PatchLevel<DIM> > level = hierarchy.getPatchLevel(ln);
    deallocatePatchData( *level );
  }
  return;
}





/*
  Deallocate patch data allocated by this class.
*/


template<int DIM>
void SinusoidalFrontTagger<DIM>::deallocatePatchData(
   hier::PatchLevel<DIM> &level )
{
  level.deallocatePatchData(d_dist_id);
  level.deallocatePatchData(d_tag_id);
  return;
}




/*
  Deallocate patch data allocated by this class.
*/
template<int DIM>
void SinusoidalFrontTagger<DIM>::computeHierarchyData(
   hier::PatchHierarchy<DIM> &hierarchy,
   double time )
{
   d_time = time;
   if ( ! d_allocate_data ) return;
   for ( int ln=0; ln<hierarchy.getNumberLevels(); ++ln ) {
      computeLevelData( hierarchy, ln, time, d_dist_id, d_tag_id );
   }
  return;
}





/*
  Compute the solution data for a level.
  Can copy data from old level (if any) to support
  initializeLevelData().
*/


template<int DIM>
void SinusoidalFrontTagger<DIM>::computeLevelData(
  const hier::PatchHierarchy<DIM> &hierarchy,
  const int ln,
  const double time,
  const int dist_id,
  const int tag_id,
  const tbox::Pointer<hier::PatchLevel<DIM> > &old_level ) const
{

  const tbox::Pointer<hier::PatchLevel<DIM> > level =
     hierarchy.getPatchLevel(ln);

  /*
    Initialize data in all patches in the level.
  */
  typename hier::PatchLevel<DIM>::Iterator pi;
  for ( pi.initialize(level); pi; pi++ ) {
    int pn = *pi;
    hier::Patch<DIM> &patch = *(level->getPatch(pn));
    tbox::Pointer<pdat::NodeData<DIM,double> > dist_data = ( dist_id >= 0 ) ?
      patch.getPatchData(dist_id) : tbox::Pointer<hier::PatchData<DIM> >(NULL);
    tbox::Pointer<pdat::CellData<DIM,int> > tag_data = ( tag_id >= 0 ) ?
      patch.getPatchData(tag_id) : tbox::Pointer<hier::PatchData<DIM> >(NULL);
    computePatchData( patch, time,
		      dist_data.getPointer(),
		      tag_data.getPointer());
  }

  return;
}







template<int DIM>
void SinusoidalFrontTagger<DIM>::overlayTagData(
   const pdat::CellData<DIM,int> &src,
   pdat::CellData<DIM,int> &dst ) const
{
   hier::Box<DIM> ibox( src.getGhostBox() * dst.getGhostBox() );
   typename hier::Box<DIM>::Iterator bi(ibox);
   for ( ; bi; bi++ ) {
      if ( src(*bi) != 0 ) {
         dst(*bi) = src(*bi);
      }
   }
   return;
}





/*
  Compute the solution data for a patch.
*/


template<int DIM>
void SinusoidalFrontTagger<DIM>::computePatchData(
  const hier::Patch<DIM> &patch,
  const double time,
  pdat::NodeData<DIM,double> *dist_data,
  pdat::CellData<DIM,int> *tag_data) const
{

   hier::Box<DIM> pbox = patch.getBox();
   tbox::Pointer<geom::CartesianPatchGeometry<DIM> > patch_geom
      = patch.getPatchGeometry();

   const double *xlo = patch_geom->getXLower();
   const double *dx = patch_geom->getDx();

   /*
     We need at least d_adaption_buffer ghost cells to compute
     the tags, but the data does not have as many ghost cells.
     So we create temporary patch data with the required "ghost"
     buffer for computing tag values.  (We could give the real
     data the required ghost cells, but that may affect the
     regridding algorithm I'm testing.)
   */
   hier::IntVector<DIM> required_tmp_buffer(d_adaption_buffer);
   pdat::NodeData<DIM,double> tmp_dist( pbox, 1, required_tmp_buffer );
   pdat::CellData<DIM,int> tmp_tag( pbox, 1, required_tmp_buffer );


   /*
     Determine what x-node-index contains the sinusoidal front.
   */

   const double wave_number = 2*3.141592654/d_period;
   double wave_offset[DIM];
   for ( int idim=0; idim<DIM; ++idim )
      wave_offset[idim] = 2*3.141592654*0;

   int i;

   hier::Box<DIM> front_box = pdat::NodeGeometry<DIM>::toNodeBox(tmp_dist.getGhostBox());
   front_box.upper(0) = front_box.lower(0);
   pdat::ArrayData<DIM,double> front_x_(front_box,1);
   MDA_Access<double,NDIM,MDA_OrderColMajor<NDIM> > front_x(
      front_x_.getPointer(0) ,
      (const int*)front_x_.getBox().lower() ,
      (const int*)front_x_.getBox().upper() );

   i = front_box.lower(0);
#if NDIM == 2
   for ( int j=front_box.lower(1); j<=front_box.upper(1); ++j ) {
      double y = xlo[1] + dx[1]*(j-pbox.lower(1));
      double siny = sin(wave_number*(y + d_init_disp[1] - d_velocity[1]*time));
      double fx = d_amplitude*siny + d_init_disp[0] + d_velocity[0]*time;
      front_x(i,j) = fx;
      // cout << i << '\t' << j << '\t' << y << '\t' << front_x(i,j) << endl;
   }
#elif NDIM == 3
   for ( int k=front_box.lower(2); k<=front_box.upper(2); ++k ) {
      double z = xlo[2] + dx[2]*(k-pbox.lower(2));
      double sinz = sin(wave_number*(z + d_init_disp[2] - d_velocity[2]*time));
      for ( int j=front_box.lower(1); j<=front_box.upper(1); ++j ) {
	 double y = xlo[1] + dx[1]*(j-pbox.lower(1));
	 double siny = sin(wave_number*(y + d_init_disp[1] + d_velocity[1]*time));
	 double fx = d_amplitude*siny*sinz + d_init_disp[0] + d_velocity[0]*time;
	 front_x(i,j,k) = fx;
	 // cout << i << '\t' << j << '\t' << k << '\t' << y << '\t' << z << '\t' << front_x(i,j,k) << endl;
      }
   }
#endif

   typename pdat::NodeData<DIM,double>::Iterator ni(tmp_dist.getGhostBox());
   for ( ; ni; ni++ ) {
      const pdat::NodeIndex<DIM> &index = *ni;
#if NDIM == 2
      tmp_dist(index) = xlo[0] + (index(0)-pbox.lower(0))*dx[0]
	 - front_x( front_box.lower(0), index(1) );
#elif NDIM == 3
      tmp_dist(index) = xlo[0] + (index(0)-pbox.lower(0))*dx[0]
	 - front_x( front_box.lower(0), index(1), index(2) );
#endif
   }
   // tmp_dist.print(tmp_dist.getBox(),0,plog);

   tmp_tag.fill(0);

   const hier::IntVector<DIM> tag_growth(d_adaption_buffer);

   typename pdat::CellData<DIM,double>::Iterator ci(tmp_dist.getGhostBox());
   for ( ; ci; ci++ ) {
      const pdat::CellIndex<DIM> &index = *ci;

      double node_x_lo = xlo[0] + (index(0)  -pbox.lower(0))*dx[0];
      double node_x_up = xlo[0] + (index(0)+1-pbox.lower(0))*dx[0];

#if NDIM == 2
      bool on_lo_side =
	 node_x_lo <= front_x(front_box.lower(0), index(1)  ) ||
	 node_x_lo <= front_x(front_box.lower(0), index(1)+1) ;
      bool on_up_side =
	 node_x_up >= front_x(front_box.lower(0), index(1)  ) ||
	 node_x_up >= front_x(front_box.lower(0), index(1)+1) ;
#elif NDIM == 3
      bool on_lo_side =
	 node_x_lo <= front_x(front_box.lower(0), index(1)  , index(2)  ) ||
	 node_x_lo <= front_x(front_box.lower(0), index(1)+1, index(2)  ) ||
	 node_x_lo <= front_x(front_box.lower(0), index(1)  , index(2)+1) ||
	 node_x_lo <= front_x(front_box.lower(0), index(1)+1, index(2)+1) ;
      bool on_up_side =
	 node_x_up >= front_x(front_box.lower(0), index(1)  , index(2)  ) ||
	 node_x_up >= front_x(front_box.lower(0), index(1)+1, index(2)  ) ||
	 node_x_up >= front_x(front_box.lower(0), index(1)  , index(2)+1) ||
	 node_x_up >= front_x(front_box.lower(0), index(1)+1, index(2)+1) ;
#endif

      if ( on_lo_side && on_up_side ) {
	 // cout << "Tagging " << index << endl;
	 hier::Box<DIM> tag_box( index, index );
	 tag_box.grow(tag_growth);
	 tag_box = tag_box * tmp_tag.getGhostBox();
	 tmp_tag.fill( 1, tag_box );
      }

   }

   /*
     Copy computed data to output.  Recall that the convention is
     to send in a NULL pointer to indicate that data is not wanted.
   */
   if ( dist_data != NULL ) {
      dist_data->copy(tmp_dist);
   }
   if ( tag_data != NULL ) {
      tag_data->copy(tmp_tag);
   }

  return;
}






template<int DIM>
int SinusoidalFrontTagger<DIM>::registerVariablesWithPlotter(
  appu::CartesianVizamraiDataWriter<DIM> &writer )
{
  /*
    Register variables with plotter.
  */
  if ( d_allocate_data ) {
    writer.registerPlotScalar("Distance to front", d_dist_id);
    writer.registerPlotScalar("Tag value", d_tag_id);
  }
  else {
    writer.registerDerivedPlotScalar("Distance to front", this);
    writer.registerDerivedPlotScalar("Tag value", this);
  }
  return 0;
}






template<int DIM>
int SinusoidalFrontTagger<DIM>::registerVariablesWithPlotter(
  appu::VisItDataWriter<DIM> &writer )
{
  /*
    Register variables with plotter.
  */
  if ( d_allocate_data ) {
    writer.registerPlotQuantity("Distance to front", "SCALAR", d_dist_id);
    writer.registerPlotQuantity("Tag value", "SCALAR", d_tag_id);
  }
  else {
    writer.registerDerivedPlotQuantity("Distance to front", "SCALAR", this,
				       1.0,
				       "NODE");
    writer.registerDerivedPlotQuantity("Tag value", "SCALAR", this);
  }
  return 0;
}






template<int DIM>
bool SinusoidalFrontTagger<DIM>::packDerivedDataIntoDoubleBuffer(
  double *buffer,
  const hier::Patch<DIM> &patch,
  const hier::Box<DIM> &region,
  const string &variable_name,
  int depth_index) const
{
   TBOX_ASSERT( d_allocate_data == false );
  if ( variable_name == "Distance to front" ) {
    pdat::NodeData<DIM,double> dist_data( patch.getBox(), 1, hier::IntVector<DIM>(0) );
    computePatchData( patch, d_time, &dist_data, NULL );
    for ( typename pdat::NodeData<DIM,double>::Iterator ci(patch.getBox()); ci; ci++ ) {
      *(buffer++) = dist_data(*ci);
    }
  }
  else if ( variable_name == "Tag value" ) {
    pdat::CellData<DIM,int> tag_data( patch.getBox() , 1, hier::IntVector<DIM>(0));
    computePatchData( patch, d_time, NULL, &tag_data );
    for ( typename pdat::CellData<DIM,double>::Iterator ci(patch.getBox()); ci; ci++ ) {
      *(buffer++) = tag_data(*ci);
    }
  }
  else {
    TBOX_ERROR("Unrecognized name " << variable_name);
  }
return true;
}


#ifdef NDIM
template class SinusoidalFrontTagger<NDIM>;
#endif
