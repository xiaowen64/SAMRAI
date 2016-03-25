//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/operators/constant/outernode/OuternodeDoubleConstantCoarsen.h $
// Package:	SAMRAI patchdata
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Constant averaging operator for node-centered double data on 
//              a  mesh.
//

#ifndef included_pdat_OuternodeDoubleConstantCoarsen
#define included_pdat_OuternodeDoubleConstantCoarsen

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_xfer_CoarsenOperator
#include "CoarsenOperator.h"
#endif


namespace SAMRAI {
   namespace pdat {

/*!
 * @brief
 * Class OuternodeDoubleConstantCoarsen implements constant
 * averaging (i.e., injection) for outernode-centered double patch data defined 
 * over a  mesh.
 *
 * It is derived from the xfer::CoarsenOperator base class.
 * The numerical operations for the averaging use FORTRAN numerical 
 * routines.
 *
 * The findCoarsenOperator() operator function returns true if the input 
 * variable is outernode-centered double, and the string is "CONSTANT_COARSEN".
 *
 * @see xfer::CoarsenOperator
 */

template<int DIM> class OuternodeDoubleConstantCoarsen
: public xfer::CoarsenOperator<DIM>
{
public:
   /*!
    * @brief Uninteresting default constructor.
    */
   OuternodeDoubleConstantCoarsen();

   /*!
    * @brief Uninteresting virtual destructor.
    */
   virtual ~OuternodeDoubleConstantCoarsen<DIM>();

   /*!
    * @brief Determine if object is for coarsening the specified
    * hier_Variable type using the given descriptive operator name.
    *
    * @return True if the variable and name string match the outernode-centered 
    * constant averaging; otherwise, return false.
    */
   bool findCoarsenOperator(const tbox::Pointer<hier::Variable<DIM> >& var,
                            const std::string &op_name) const; 

   /*!
    * @brief Return descriptive name string identifier of this
    * coarsening operator.
    *
    * @return descriptive name string identifier of this coarsening operator.
    */
   const std::string& getOperatorName() const;

   /*!
    * @brief Give the operator priority.
    *
    * The priority of outernode-centered constant averaging is 0.
    * It will be performed before any user-defined coarsen operations. 
    */
   int getOperatorPriority() const;

   /*!
    * @brief Give the operator stencil width.
    *
    * The stencil width of the constant averaging operator is the vector of 
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector<DIM> getStencilWidth() const;

   /*!
    * @brief Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the outernode-centered double constant 
    * averaging operator.
    *
    * Coarsening is performed on the intersection of 
    * the destination patch and the coarse box.  It is assumed that the 
    * fine patch contains sufficient data for the stencil width of the 
    * coarsening operator.
    */
   void coarsen(hier::Patch<DIM>& coarse,
                const hier::Patch<DIM>& fine,
                const int dst_component,
                const int src_component,
                const hier::Box<DIM>& coarse_box,
                const hier::IntVector<DIM>& ratio) const;

private:
   std::string d_name_id;

};

}
}

#endif
