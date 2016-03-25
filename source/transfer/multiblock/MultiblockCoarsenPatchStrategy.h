//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/transfer/multiblock/MultiblockCoarsenPatchStrategy.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Strategy interface to user routines for coarsening AMR data.
//
 
#ifndef included_xfer_MultiblockCoarsenPatchStrategy
#define included_xfer_MultiblockCoarsenPatchStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_xfer_CoarsenPatchStrategy
#include "CoarsenPatchStrategy.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

#ifndef MULTIBLOCK_UNDEFINED_BLOCK_NUMBER
#define MULTIBLOCK_UNDEFINED_BLOCK_NUMBER (-1)
#endif

namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Class MultiblockCoarsenPatchStrategy is an abstract base class that
 * provides interfaces between coarsen operations in the multiblock 
 * classes and user's problem specific routines.  In general, it mimics
 * the standard CoarsenPatchStrategy class but implements some additional
 * methods that supply information about the blocks.
 *
 * @see xfer::CoarsenPatchStrategy
 * @see xfer::MultiblockCoarsenSchedule
 */
 
template<int DIM> class MultiblockCoarsenPatchStrategy
:
public virtual tbox::DescribedClass,
public virtual xfer::CoarsenPatchStrategy<DIM>
{
public:
   /*!
    * The constructor for coarsen strategy does nothing interesting.
    */
   MultiblockCoarsenPatchStrategy();
 
   /*!
    * The virtual destructor for coarsen strategy does nothing interesting.
    */
   virtual ~MultiblockCoarsenPatchStrategy<DIM>();

   /*!
    * Return maximum stencil width needed over all user-defined
    * data coarsening operations.  This is needed to
    * determine the correct coarsening data dependencies.
    */
   virtual hier::IntVector<DIM> getMultiblockCoarsenOpStencilWidth() = 0;

   /*!
    * Perform user-defined coarsening operations.  This member function
    * is called before standard coarsening operations (expressed using
    * concrete subclasses of the MultiblockCoarsenOperator<DIM> base class).  
    * The preprocess function should move data from the source components 
    * on the fine patch into the source components on the coarse patch 
    * in the specified coarse box region.   Recall that the source components
    * are specified in calls to the registerCoarsen() function in the 
    * MultiblockCoarsenAlgorithm<DIM> class.
    *
    * @param coarse      Coarse patch containing destination data.
    * @param fine        Fine patch containing source data.
    * @param coarse_box  hier::Box region on coarse patch into which data is coarsened. 
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void preprocessCoarsen(
      hier::Patch<DIM>& coarse,
      const hier::Patch<DIM>& fine,
      const hier::Box<DIM>& coarse_box,
      const hier::IntVector<DIM>& ratio)
      {
         NULL_USE(coarse);
         NULL_USE(fine);
         NULL_USE(coarse_box);
         NULL_USE(ratio);
         return;
      }

   /*!
    * Perform user-defined coarsening operations.  This member function
    * is called after standard coarsening operations (expressed using
    * concrete subclasses of the CoarsenOperator<DIM> base class).  
    * The postprocess function should move data from the source components on 
    * the fine patch into the source components on the coarse patch in the
    * specified coarse box region.  Recall that the source components are 
    * specified in calls to the registerMultiblockCoarsen() function in the 
    * MultiblockCoarsenAlgorithm<DIM> class.
    *
    * @param coarse      Coarse patch containing destination data.
    * @param fine        Fine patch containing source data.
    * @param coarse_box  hier::Box region on coarse patch into which data is copied.
    * @param ratio       Integer vector containing ratio
    */
   virtual void postprocessCoarsen(
      hier::Patch<DIM>& coarse,
      const hier::Patch<DIM>& fine,
      const hier::Box<DIM>& coarse_box,
      const hier::IntVector<DIM>& ratio)
      {
         NULL_USE(coarse);
         NULL_USE(fine);
         NULL_USE(coarse_box);
         NULL_USE(ratio);
         return;
      }

   /*
    * Set the multiblock block number.
    */
   virtual void setCoarsenBlockNumber(const int block_number)
   {
      d_block_number = block_number;
   }

   /*
    * Get the multiblock block number.
    */
   virtual int getCoarsenBlockNumber()
   {
      return(d_block_number);
   }

   /*
    * Clear the multiblock block number.
    */
   virtual void clearCoarsenBlockNumber()
   {
      d_block_number = MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;
   }

protected:

   int  d_block_number;

private:

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "MultiblockCoarsenPatchStrategy.C"
#endif
