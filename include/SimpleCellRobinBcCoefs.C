#ifndef included_solv_SimpleCellRobinBcCoefs_C
#define included_solv_SimpleCellRobinBcCoefs_C

/*
 * File:         $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/solvers/poisson/SimpleCellRobinBcCoefs.C $
 * Package:      SAMRAI solvers
 * Copyright:    (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:     $LastChangedRevision: 1846 $
 * Modified:     $LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
 * Description:  Level solver for diffusion-like elliptic problems.
 */

#include "CartesianPatchGeometry.h"
#include "ArrayDataBasicOps.h"
#include "CellData.h"
#include "OuterfaceData.h"
#include "SideData.h"
#include "SimpleCellRobinBcCoefs.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"


#ifndef NULL
#define NULL (0)
#endif

#define DIRICHLET 0
#define NEUMANN 1
#define MIXED 2


namespace SAMRAI {
    namespace solv {

/*
*************************************************************************
*                                                                       *
* Construct an unitialized boundary specification.                      *
*                                                                       *
*************************************************************************
*/

template<int DIM>  SimpleCellRobinBcCoefs<DIM>::SimpleCellRobinBcCoefs(
   const std::string& object_name)
   : d_object_name(object_name),
     d_hierarchy(NULL),
     d_ln_min(-1),
     d_ln_max(-1),
     d_flux_id(-1),
     d_flag_id(-1),
     d_dirichlet_data_id(-1),
     d_diffusion_coef_id(-1)
{
   t_set_bc_coefs = tbox::TimerManager::getManager()->
      getTimer("solv::SimpleCellRobinBcCoefs::setBcCoefs()");

   return;
}




template<int DIM>  SimpleCellRobinBcCoefs<DIM>::~SimpleCellRobinBcCoefs()
{
  return;
}



template<int DIM> void SimpleCellRobinBcCoefs<DIM>::setHierarchy(
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int ln_min,
   const int ln_max)
{
   d_hierarchy = hierarchy;
   d_ln_min = ln_min;
   d_ln_max = ln_max;
   d_dirichlet_data.setNull();
   d_dirichlet_data_pos.setNull();
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed.\n");
   }
#endif
   if ( d_ln_min == -1 ) {
      d_ln_min = 0;
   }
   if ( d_ln_max == -1 ) {
      d_ln_min = d_hierarchy->getFinestLevelNumber();
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max ) {
      TBOX_ERROR(d_object_name
                 << ": Bad range of levels in setHierarchy().\n");
   }
#endif
   return;
}



template<int DIM> void SimpleCellRobinBcCoefs<DIM>::setBoundaries(
   const std::string &boundary_type,
   const int fluxes,
   const int flags,
   int *bdry_types)
{

   int k;

   if (boundary_type == "Dirichlet") {
      d_flux_id = -1;
      d_flag_id = -1;
      for (k=0;k<2*DIM;k++) {
         d_bdry_types[k] = DIRICHLET;
      }
   } else if (boundary_type == "Neumann") {
      for (k=0;k<2*DIM;k++) {
         d_bdry_types[k] = NEUMANN;
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      if (fluxes < 0) {
         TBOX_ERROR(d_object_name << ": bad flux patch data index ("
                    << fluxes << ") for Neumann boundary condition.\n");
      }
#endif
      d_flux_id = fluxes;
      d_flag_id = -1;
   } else if (boundary_type == "Mixed") {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (fluxes < 0) {
         TBOX_ERROR(d_object_name << ": bad flux patch data index ("
                    << fluxes << ") for Mixed boundary condition.\n");
      }
      if (flags  < 0) {
         TBOX_ERROR(d_object_name << ": bad flag patch data index ("
                    << flags << ") for Mixed boundary condition.\n");
      }
#endif
      d_flux_id = fluxes;
      d_flag_id = flags;
      if (bdry_types != NULL) {
         for (k=0;k<2*DIM;k++) {
            d_bdry_types[k] = bdry_types[k];
         }
      } else {
         for (k=0;k<2*DIM;k++) {
            d_bdry_types[k] = MIXED;
         }
       } 
   } else {
      TBOX_ERROR(d_object_name << ": Non-existing case of\n"
                 << "boundary_type in PoissonSolver<DIM>::setBoundaries()");
   }

}





/*
************************************************************************
*                                                                      *
* Set the bc coefficients according to information received in the     *
* call to setBoundaries.                                               *
*                                                                      *
* Do a lot of error checking before hand.                              *
*                                                                      *
* For Dirichlet, we need the Dirichlet boundary                        *
* values stored in the ghost cells of the solution.                    *
*                                                                      *
* For Neumann bc, we need the flux data and the                        *
* diffusion coefficient to determine the required                      *
* normal gradient of the solution.  However, the                       *
* diffusion coefficient is assumed to be 1.0 if                        *
* left unspecified.                                                    *
*                                                                      *
* For mixed bc, we need the flag stating whether                       *
* Dirichlet or Neumann at any face, in addition to                     *
* Dirichlet and Neumann data.                                          *
*                                                                      *
************************************************************************
*/

template<int DIM> void SimpleCellRobinBcCoefs<DIM>::setBcCoefs (
   tbox::Pointer<pdat::ArrayData<DIM,double> > &acoef_data ,
   tbox::Pointer<pdat::ArrayData<DIM,double> > &bcoef_data ,
   tbox::Pointer<pdat::ArrayData<DIM,double> > &gcoef_data ,
   const tbox::Pointer< hier::Variable<DIM> > &variable ,
   const hier::Patch<DIM> &patch ,
   const hier::BoundaryBox<DIM> &bdry_box ) const
{
   NULL_USE(variable);

   t_set_bc_coefs->start();

   const int ln = patch.getPatchLevelNumber();
   const int pn = patch.getPatchNumber();
   const int location_index = bdry_box.getLocationIndex();

   tbox::Pointer< hier::PatchData<DIM> > patch_data;
   tbox::Pointer<pdat::OuterfaceData<DIM,double> > flux_data_ptr;
   tbox::Pointer<pdat::SideData<DIM,double> > diffcoef_data_ptr;
   tbox::Pointer<pdat::OuterfaceData<DIM,int> > flag_data_ptr;

#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !gcoef_data.isNull() ) {
      if ( d_bdry_types[location_index] == DIRICHLET 
           || d_bdry_types[location_index] == MIXED ) {
         /*
          * For Dirichlet and mixed boundaries, we use cached data
          * to get the Dirichlet value.  Data specific to the
          * d_hierarchy has been cached by cacheDirichletData().
          * This implementation can only set Dirichlet coefficients
          * when the patch is in the correct hierarchy.
          */
         if ( ! patch.inHierarchy() ) {
            TBOX_ERROR(d_object_name << ": patch is not in any hierarchy.\n"
                       << "SimpleCellRobinBcCoefs<DIM> can only set\n"
                       << "boundary coefficients for patches in\n"
                       << "the same hierarchy as cached\n"
                       << "Dirichlet coefficients.");
         }
         tbox::Pointer<hier::PatchLevel<DIM> > level = d_hierarchy->getPatchLevel(ln);
         if ( level->getPatch(pn)->getBox() != patch.getBox() ) {
            TBOX_ERROR(d_object_name << ": patch is not in the hierarchy\n"
                       << "of cached boundary data.\n"
                       << "SimpleCellRobinBcCoefs<DIM> can only set\n"
                       << "boundary coefficients for patches in\n"
                       << "the same hierarchy as cached\n"
                       << "Dirichlet coefficients.");
         }
      }
      if ( d_bdry_types[location_index] == NEUMANN 
           || d_bdry_types[location_index] == MIXED ) {
         patch_data = patch.getPatchData( d_flux_id );
         flux_data_ptr = patch_data;
         if ( !patch_data ) {
            TBOX_ERROR(d_object_name << ": Flux data (patch data id = "
                       << d_flux_id << ") does not exist.");
         }
         if ( !flux_data_ptr ) {
            TBOX_ERROR(d_object_name << ": Flux data (patch data id = "
                       << d_flux_id << ") is not outerface double data.");
         }
         if ( d_diffusion_coef_id != -1 ) {
            patch_data = patch.getPatchData( d_diffusion_coef_id );
            diffcoef_data_ptr = patch_data;
            if ( !patch_data ) {
               TBOX_ERROR(d_object_name << ": Diffusion coefficient data\n"
                          "(patch data id = " << d_diffusion_coef_id
                          << ") does not exist.");
            }
            if ( !diffcoef_data_ptr ) {
               TBOX_ERROR(d_object_name << ": Diffusion coefficient data\n"
                          "(patch data id = " << d_diffusion_coef_id
                          << ") is not side double data.");
            }
         }
      }
   }
   if ( !acoef_data.isNull() ) {
      if ( d_bdry_types[location_index] == MIXED ) {
         patch_data = patch.getPatchData( d_flag_id );
         flag_data_ptr = patch.getPatchData( d_flag_id );
         if ( !patch_data ) {
            TBOX_ERROR(d_object_name << ": Flags data (patch data id = "
                       << d_flag_id << ") does not exist.");
         }
         if ( !flag_data_ptr ) {
            TBOX_ERROR(d_object_name << ": Flags data (patch data id = "
                       << d_flag_id << ") is not outerface int data.");
         }
      }
   }
#endif


   int bn;


   if ( d_bdry_types[location_index] == DIRICHLET ) {

      if ( !acoef_data.isNull() ) {
         acoef_data->fill( 1.0 );
      }
      if ( !bcoef_data.isNull() ) {
         bcoef_data->fill( 0.0 );
      }

      if ( !gcoef_data.isNull() ) {

         tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
            patch.getPatchGeometry();
         const tbox::Array< hier::BoundaryBox<DIM> > &codim1_boxes =
            pg->getCodimensionBoundary(1);
         /*
	   Search for cached boundary box containing current boundary box.
          */
         for ( bn=0; bn<codim1_boxes.getSize(); ++bn ) {
            const hier::BoundaryBox<DIM> &cdb = codim1_boxes[bn];
            if ( bdry_box.getLocationIndex() == cdb.getLocationIndex()
              && bdry_box.getBox().lower() >= cdb.getBox().lower()
              && bdry_box.getBox().upper() <= cdb.getBox().upper()
                 ) break;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         if ( bn == codim1_boxes.getSize() ) {
            TBOX_ERROR(d_object_name << " cannot find cached Dirichlet data.\n"
                       << "This is most likely caused by not calling\n"
                       << "SimpleCellRobinBcCoefs<DIM>::cacheDirichletData()\n"
                       << "after the hierarchy changed.\n");
         }
#endif
         int position = d_dirichlet_data_pos[ln][pn] + bn;
         gcoef_data->copy( *d_dirichlet_data[position],
                           d_dirichlet_data[position]->getBox(),
                           hier::IntVector<DIM>(0) );

      }
   }

   else if ( d_bdry_types[location_index] == NEUMANN ) {

      if ( !acoef_data.isNull() ) {
         acoef_data->fill( 0.0 );
      }
      if ( !bcoef_data.isNull() ) {
         bcoef_data->fill( 1.0 );
      }

      if ( !gcoef_data.isNull() ) {
         flux_data_ptr = patch.getPatchData(d_flux_id);
	 pdat::OuterfaceData<DIM, double> &flux_data(*flux_data_ptr);
	 const int axis = location_index/2;
	 const int face = location_index%2;
         pdat::ArrayData<DIM,double> &g = *gcoef_data;
         pdat::ArrayDataIterator<DIM> ai(g.getBox());
	 hier::Index<DIM> offset_to_inside(0);
	 if ( face != 0 ) offset_to_inside(axis) = -1;
         if ( d_diffusion_coef_id == -1 ) {
            for ( ; ai; ai++ ) {
	       pdat::FaceIndex<DIM> fi( ai()+offset_to_inside, axis, face );
               g(*ai,0) = flux_data(fi,face)/d_diffusion_coef_constant;
	       tbox::plog << location_index << '\t' << g(*ai,0) << '\n';
            }
         }
         else {
            diffcoef_data_ptr = patch.getPatchData(d_diffusion_coef_id);
            const pdat::ArrayData<DIM,double> &diffcoef_array_data
               = diffcoef_data_ptr->getArrayData(axis);
            for ( ; ai; ai++ ) {
	       pdat::FaceIndex<DIM> fi( ai()+offset_to_inside, axis, face );
               g(*ai,0) = flux_data(fi,face)/diffcoef_array_data(*ai,0);
	       tbox::plog << location_index << '\t' << g(*ai,0) << '\n';
            }
         }
      }

   }

   else if ( d_bdry_types[location_index] == MIXED ) {

      const int axis = location_index/2;
      const int face = location_index%2;
      flag_data_ptr = patch.getPatchData(d_flag_id); 
      pdat::OuterfaceData<DIM, int> &flag_data(*flag_data_ptr);
      hier::Index<DIM> offset_to_inside(0);
      if ( face != 0 ) offset_to_inside(axis) = -1;

      if ( !acoef_data.isNull() ) {
         pdat::ArrayData<DIM,double> &a = *acoef_data;
         pdat::ArrayDataIterator<DIM> ai(a.getBox());
         for ( ; ai; ai++ ) {
	    pdat::FaceIndex<DIM> fi( ai()+offset_to_inside, axis, face );
            if ( flag_data(fi,face) == 0 ) {
               a(*ai,0) = 1.0;
            }
            else {
               a(*ai,0) = 0.0;
            }
         }
      }

      if ( !bcoef_data.isNull() ) {
         pdat::ArrayData<DIM,double> &b = *bcoef_data;
         pdat::ArrayDataIterator<DIM> bi(b.getBox());
         for ( ; bi; bi++ ) {
	    pdat::FaceIndex<DIM> fi( bi()+offset_to_inside, axis, face );
            if ( flag_data(fi,face) == 0 ) {
               b(*bi,0) = 0.0;
            }
            else {
               b(*bi,0) = 1.0;
            }
         }
      }

      if ( !gcoef_data.isNull() ) {
         tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
            patch.getPatchGeometry();
         const tbox::Array< hier::BoundaryBox<DIM> > &codim1_boxes =
            pg->getCodimensionBoundary(1);
         /*
	   Search for cached boundary box containing current boundary box.
          */
         for ( bn=0; bn<codim1_boxes.getSize(); ++bn ) {
            const hier::BoundaryBox<DIM> &cdb = codim1_boxes[bn];
            if ( bdry_box.getLocationIndex() == cdb.getLocationIndex()
              && bdry_box.getBox().lower() >= cdb.getBox().lower()
              && bdry_box.getBox().upper() <= cdb.getBox().upper()
                 ) break;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         if ( bn == codim1_boxes.getSize() ) {
            TBOX_ERROR(d_object_name << " cannot find cached Dirichlet data.\n"
                       << "This is most likely caused by not calling\n"
                       << "SimpleCellRobinBcCoefs<DIM>::cacheDirichletData() after the\n"
                       << "hierarchy changed.\n");
         }
#endif
         int position = d_dirichlet_data_pos[ln][pn] + bn;
         const pdat::ArrayData<DIM,double> &dirichlet_array_data =
            *d_dirichlet_data[position];
         diffcoef_data_ptr = patch.getPatchData(d_diffusion_coef_id);
         pdat::ArrayData<DIM, double> &g = *gcoef_data;
	 pdat::OuterfaceData<DIM, double> &flux_data(*flux_data_ptr);
         pdat::ArrayDataIterator<DIM> ai(g.getBox());
         for ( ; ai; ai++ ) {
	    pdat::FaceIndex<DIM> fi( ai()+offset_to_inside, axis, face );
            if ( flag_data(fi,face) == 0 ) {
               g(*ai,0) = dirichlet_array_data(*ai,0);
            }
            else {
	       pdat::FaceIndex<DIM> fi2( ai()+offset_to_inside, axis, face );
               if ( d_diffusion_coef_id == -1 ) {
                  g(*ai,0) = flux_data(fi2,face)/d_diffusion_coef_constant;
               }
               else {
                  g(*ai,0) = flux_data(fi2,face)/
		    diffcoef_data_ptr->getArrayData(axis)(*ai,0);
               }
            }
         }
      }

   }

   t_set_bc_coefs->stop();

   return;
}



/*
***********************************************************************
* This class cannot set coefficients for boundary boxes that extend   *
* past the patch in the direction parallel to the boundary,           *
* because it relies on data, such as pdat::OutersideData<DIM>,             *
* that does not extend.                                               *
***********************************************************************
*/
template<int DIM> hier::IntVector<DIM> SimpleCellRobinBcCoefs<DIM>::numberOfExtensionsFillable() const
{
   return hier::IntVector<DIM>(0);
}



/*
************************************************************************
*                                                                      *
* Copy and save cell-centered Dirichlet data in ghost cells.           *
* For each boundary box in the hierarchy, we create a pdat::ArrayData<DIM>  *
* object indexed on the side indices corresponding to boundary boxes.  *
* The ghost-cell-centered data is shifted to the side indices and      *
* saved in the pdat::ArrayData<DIM> objects.                                *
*                                                                      *
* First, loop through the hierarchy to compute how many                *
* pdat::ArrayData<DIM> objects we need and the position of each one.        *
*                                                                      *
* Second, allocate the pdat::ArrayData<DIM> objects.                        *
*                                                                      *
* Third, loop through the hierarchy again to allocate the data in each *
* pdat::ArrayData<DIM> object and cache the ghost data.                     *
*                                                                      *
* The position of the appropriate boundary box bn of patch pn          *
* of level ln is d_dirichlet_data_pos[ln][pn]+bn                       *
*                                                                      *
************************************************************************
*/
template<int DIM> void SimpleCellRobinBcCoefs<DIM>::cacheDirichletData( int dirichlet_data_id )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( !d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": hierarchy has not been set.\n"
                 << "use setHierarchy() to set the hierarchy before\n"
                 << "caching boundary ghost cell data.\n");
   }
#endif
   d_dirichlet_data.setNull();
   d_dirichlet_data_pos.setNull();
   int i, ln, pn, bn, position, n_reqd_boxes=0;
   d_dirichlet_data_pos.resizeArray(d_ln_max+1);
   for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
      hier::PatchLevel<DIM> &level = (hier::PatchLevel<DIM>&)
	 *d_hierarchy->getPatchLevel(ln);
      const int num_patches = level.getNumberOfPatches();
      d_dirichlet_data_pos[ln].resizeArray(num_patches);
      for ( i=0; i<num_patches; ++i ) {
         d_dirichlet_data_pos[ln][i] = -1;
      }
      typename hier::PatchLevel<DIM>::Iterator pi(level);
      for ( ; pi; pi++ ) {
         pn = *pi;
         hier::Patch<DIM> &patch = *level.getPatch(pn);
         tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
            patch.getPatchGeometry();
         const tbox::Array< hier::BoundaryBox<DIM> > &codim1_boxes =
            pg->getCodimensionBoundary(1);
         d_dirichlet_data_pos[ln][pn] = n_reqd_boxes;
         n_reqd_boxes += codim1_boxes.getSize();
      }
   }
   d_dirichlet_data.resizeArray(n_reqd_boxes);
   for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
      hier::PatchLevel<DIM> &level = (hier::PatchLevel<DIM>&)
	 *d_hierarchy->getPatchLevel(ln);
      typename hier::PatchLevel<DIM>::Iterator pi(level);
      for ( ; pi; pi++ ) {
         pn = *pi;
         hier::Patch<DIM> &patch = *level.getPatch(pn);
         tbox::Pointer<pdat::CellData<DIM,double> > cell_data =
            patch.getPatchData(dirichlet_data_id);
         tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
            patch.getPatchGeometry();
         const tbox::Array< hier::BoundaryBox<DIM> > &codim1_boxes =
            pg->getCodimensionBoundary(1);
         for ( bn=0; bn<codim1_boxes.getSize(); ++bn ) {
            const hier::BoundaryBox<DIM> &bdry_box = codim1_boxes[bn];
            position = d_dirichlet_data_pos[ln][pn] + bn;
            hier::Box<DIM> databox = makeSideBoundaryBox(bdry_box);
            d_dirichlet_data[position] =
              new pdat::ArrayData<DIM,double>( databox,
                                           1,
                                           tbox::Pointer<tbox::Arena>(0));
            pdat::ArrayData<DIM,double> &array_data = *d_dirichlet_data[position];
            hier::IntVector<DIM> shift_amount(0);
            const int location_index = bdry_box.getLocationIndex();
            if ( location_index%2 == 0 ) shift_amount[location_index/2] = 1;
            array_data.copy(cell_data->getArrayData(),
                            databox,
                            shift_amount);
         }
      }
   }
   return;
}



/*
************************************************************************
*                                                                      *
* Reverse action of cacheDirichletData by copying cached data back     *
* into the ghost cells.                                                *
*                                                                      *
* The cached data is not deallocated.                                  *
*                                                                      *
************************************************************************
*/
template<int DIM> void SimpleCellRobinBcCoefs<DIM>::restoreDirichletData( int dirichlet_data_id )
{
   if ( d_dirichlet_data_pos.isNull() ) {
      TBOX_ERROR(d_object_name << ".restoreDirichletData(): Dirichlet\n"
                 << "data has not been set.\n");
   }
   int ln, pn, bn, position;
   for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
      hier::PatchLevel<DIM> &level = (hier::PatchLevel<DIM>&)
	 *d_hierarchy->getPatchLevel(ln);
      typename hier::PatchLevel<DIM>::Iterator pi(level);
      for ( ; pi; pi++ ) {
         pn = *pi;
         hier::Patch<DIM> &patch = *level.getPatch(pn);
         tbox::Pointer<pdat::CellData<DIM,double> > cell_data =
            patch.getPatchData(dirichlet_data_id);
         tbox::Pointer< geom::CartesianPatchGeometry<DIM> > pg =
            patch.getPatchGeometry();
         const tbox::Array< hier::BoundaryBox<DIM> > &codim1_boxes =
            pg->getCodimensionBoundary(1);
         for ( bn=0; bn<codim1_boxes.getSize(); ++bn ) {
            const hier::BoundaryBox<DIM> &bdry_box = codim1_boxes[bn];
            position = d_dirichlet_data_pos[ln][pn] + bn;
            pdat::ArrayData<DIM,double> &array_data = *d_dirichlet_data[position];
            hier::IntVector<DIM> shift_amount(0);
            const int location_index = bdry_box.getLocationIndex();
            hier::Box<DIM> dst_box=array_data.getBox();
            if ( location_index%2 == 0 ) {
               shift_amount[location_index/2] = -1;
               dst_box.shift(location_index/2,-1);
            }
            cell_data->getArrayData().copy(array_data,
                                           dst_box,
                                           shift_amount);
         }
      }
   }
   return;
}




template<int DIM> void SimpleCellRobinBcCoefs<DIM>::setDiffusionCoefId( int diffusion_coef_id )
{
  d_diffusion_coef_id = diffusion_coef_id;
  d_diffusion_coef_constant = 0.0;
  return;
}




template<int DIM> void SimpleCellRobinBcCoefs<DIM>::setDiffusionCoefConstant( double diffusion_coef_constant )
{
  d_diffusion_coef_constant = diffusion_coef_constant;
  d_diffusion_coef_id = -1;
  return;
}



/*
************************************************************************
* Make surface box on boundary using standard boundary box             *
************************************************************************
*/

template<int DIM> hier::Box<DIM> SimpleCellRobinBcCoefs<DIM>::makeSideBoundaryBox(
   const hier::BoundaryBox<DIM> &boundary_box ) const
{
   if ( boundary_box.getBoundaryType() != 1 ) {
      TBOX_ERROR(d_object_name << ": CartesianRobinBcHelper<DIM>::makeSideBoundaryBox called with\n"
                 << "improper boundary box\n"
                 << "for " << d_object_name );
   }
   hier::Box<DIM> face_indices = boundary_box.getBox();
   int location_index = boundary_box.getLocationIndex();
   if ( location_index%2 == 0 ) {
      /*
       * On the min index side, the face indices are one higher
       * than the boundary cell indices, in the direction normal
       * to the boundary.
       */
      face_indices.shift(location_index/2,1);
   }
   return face_indices;
}


}
}
#endif
