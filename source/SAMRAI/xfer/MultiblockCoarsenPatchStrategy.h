/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data. 
 *
 ************************************************************************/

#ifndef included_xfer_MultiblockCoarsenPatchStrategy
#define included_xfer_MultiblockCoarsenPatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/tbox/Utilities.h"

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

class MultiblockCoarsenPatchStrategy:
   public virtual tbox::DescribedClass,
   public xfer::CoarsenPatchStrategy
{
public:
   /*!
    * The constructor for coarsen strategy does nothing interesting.
    */
   explicit MultiblockCoarsenPatchStrategy(
      const tbox::Dimension& dim);

   /*!
    * The virtual destructor for coarsen strategy does nothing interesting.
    */
   virtual ~MultiblockCoarsenPatchStrategy();

   /*!
    * Return maximum stencil width needed over all user-defined
    * data coarsening operations.  This is needed to
    * determine the correct coarsening data dependencies.
    */
   virtual hier::IntVector
   getMultiblockCoarsenOpStencilWidth() = 0;

   /*!
    * Perform user-defined coarsening operations.  This member function
    * is called before standard coarsening operations (expressed using
    * concrete subclasses of the MultiblockCoarsenOperator<DIM> base class).
    * The preprocess function should move data from the source components
    * on the fine patch into the source components on the coarse patch
    * in the specified coarse box region.   Recall that the source components
    * are specified in calls to the registerCoarsen() function in the
    * MultiblockCoarsenAlgorithm class.
    *
    * @param coarse      Coarse patch containing destination data.
    * @param fine        Fine patch containing source data.
    * @param coarse_box  hier::Box region on coarse patch into which data is coarsened.
    * @param ratio       Integer vector containing ratio relating index space
    *                    between coarse and fine patches.
    */
   virtual void preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   /*!
    * Perform user-defined coarsening operations.  This member function
    * is called after standard coarsening operations (expressed using
    * concrete subclasses of the CoarsenOperator<DIM> base class).
    * The postprocess function should move data from the source components on
    * the fine patch into the source components on the coarse patch in the
    * specified coarse box region.  Recall that the source components are
    * specified in calls to the registerMultiblockCoarsen() function in the
    * MultiblockCoarsenAlgorithm class.
    *
    * @param coarse      Coarse patch containing destination data.
    * @param fine        Fine patch containing source data.
    * @param coarse_box  hier::Box region on coarse patch into which data is copied.
    * @param ratio       Integer vector containing ratio
    */
   virtual void postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(coarse);
      NULL_USE(fine);
      NULL_USE(coarse_box);
      NULL_USE(ratio);
   }

   /*
    * Set the multiblock block number.
    */
   virtual void setCoarsenBlockNumber(
      const int block_number)
   {
      d_block_number = block_number;
   }

   /*
    * Get the multiblock block number.
    */
   virtual int getCoarsenBlockNumber()
   {
      return d_block_number;
   }

   /*
    * Clear the multiblock block number.
    */
   virtual void clearCoarsenBlockNumber()
   {
      d_block_number = MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;
   }

protected:
   int d_block_number;

private:
   /*
    * Static integer constant describing an undefined block.
    */
   static const int MULTIBLOCK_UNDEFINED_BLOCK_NUMBER;

};

}
}
#endif
