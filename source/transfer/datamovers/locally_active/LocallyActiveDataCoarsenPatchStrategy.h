//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/locally_active/LocallyActiveDataCoarsenPatchStrategy.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Strategy interface to user routines for coarsening locally-active AMR data.
//
 
#ifndef included_xfer_LocallyActiveDataCoarsenPatchStrategy
#define included_xfer_LocallyActiveDataCoarsenPatchStrategy

#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/DescribedClass.h"
#include "tbox/List.h"

namespace SAMRAI {
   namespace xfer {

/*!
 * @brief Class LocallyActiveDataCoarsenPatchStrategy is an abstract 
 * base class that defines a Strategy pattern interface for operations that 
 * a user may supply for application-specific coarsening of locally-active data 
 * between two levels in an AMR patch hierarchy.  This interface class is similar to
 * the CoarsenPatchStrategy interface, except that it is used to treat locally-active
 * patch data.  A concrete subclass must define three member functions to perform
 * the following tasks:
 *
 * - define maximum stencil width for user-defined coarsen operations
 * - preprocess the coarsening
 * - postprocess the coarsening
 *
 * Note that the preprocess member function is called before standard data coarsening 
 * using CoarsenOperators and the postprocessor member function is called afterwards.
 *
 * @see xfer::LocallyActiveDataCoarsenAlgorithm
 * @see xfer::LocallyActiveDataCoarsenSchedule
 */

template<int DIM> class 
LocallyActiveDataCoarsenPatchStrategy
:
public virtual tbox::DescribedClass
{
public:
   /*!
    * The constructor for the coarsen strategy does nothing interesting.
    */
   LocallyActiveDataCoarsenPatchStrategy();
 
   /*!
    * The virtual destructor for the coarsen strategy does nothing interesting.
    */
   virtual ~LocallyActiveDataCoarsenPatchStrategy();

   /*!
    * Return maximum stencil width needed over all user-defined
    * data coarsening operations.  This is needed to
    * determine the correct coarsening data dependencies.
    */
   virtual hier::IntVector<DIM> getCoarsenOpStencilWidth() const = 0;

   /*!
    * Perform user-defined coarsening operations.  This member function
    * is called before standard coarsening operations (expressed using
    * concrete subclasses of the CoarsenOperator<DIM> base class).  
    * The preprocess function should move data from the source components 
    * on the fine patch into the source components on the coarse patch 
    * in the specified coarse box region.   Recall that the source components
    * are specified in calls to the registerCoarsen() function in the 
    * CoarsenAlgorithm<DIM> class.
    *
    * @param coarse        Coarse patch containing destination data.
    * @param fine          Fine patch containing source data.
    * @param src_data_ids  List of integer patch data ids to preprocess.
    * @param coarse_box    Box region on coarse patch into which data is coarsened. 
    * @param ratio         Integer vector containing ratio relating index space
    *                      between coarse and fine patches.
    */
   virtual void preprocessCoarsen(
      hier::Patch<DIM>& coarse,
      const hier::Patch<DIM>& fine,
      const tbox::List<int>& src_data_ids,
      const hier::Box<DIM>& coarse_box,
      const hier::IntVector<DIM>& ratio) = 0;

   /*!
    * Perform user-defined coarsening operations.  This member function
    * is called after standard coarsening operations (expressed using
    * concrete subclasses of the CoarsenOperator<DIM> base class).  
    * The postprocess function should move data from the source components on 
    * the fine patch into the source components on the coarse patch in the
    * specified coarse box region.  Recall that the source components are 
    * specified in calls to the registerCoarsen() function in the 
    * CoarsenOperator<DIM> class.
    *
    * @param coarse        Coarse patch containing destination data.
    * @param fine          Fine patch containing source data.
    * @param src_data_ids  List of integer patch data ids to preprocess.
    * @param coarse_box    Box region on coarse patch into which data is copied.
    * @param ratio         Integer vector containing ratio
    */
   virtual void postprocessCoarsen(
      hier::Patch<DIM>& coarse,
      const hier::Patch<DIM>& fine,
      const tbox::List<int>& src_data_ids,
      const hier::Box<DIM>& coarse_box,
      const hier::IntVector<DIM>& ratio) = 0;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataCoarsenPatchStrategy.C"
#endif
