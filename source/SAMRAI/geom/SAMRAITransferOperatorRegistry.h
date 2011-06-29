/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Singleton registry for all tranfer operators.
 *
 ************************************************************************/

#ifndef included_hier_SAMRAITransferOperatorRegistry
#define included_hier_SAMRAITransferOperatorRegistry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/TransferOperatorRegistry.h"

#include <string>

namespace SAMRAI {
namespace geom {

/*!
 * @brief Class SAMRAITransferOperatorRegistry is intended to implement the
 * buildCoarsenOperator, buildRefineOperator and buildTimeInterpolateOperator
 * member functions of TransferOperatorRegistry.

 * @see hier::TransferOperatorRegistry
 */
  class SAMRAITransferOperatorRegistry : public hier::TransferOperatorRegistry
{
public:
   /*!
    * @brief Set the state of the hier::TransferOperatorRegistry members.
    *
    * @param[in]  dim
    */
   SAMRAITransferOperatorRegistry(const tbox::Dimension& dim);

   /*!
    * @brief Destructor
    */
   ~SAMRAITransferOperatorRegistry();

protected:
   tbox::Pointer<hier::CoarsenOperator> buildCoarsenOperator(
      const tbox::Pointer<hier::Variable>& var,
      const std::string& op_name);

   tbox::Pointer<hier::RefineOperator> buildRefineOperator(
      const tbox::Pointer<hier::Variable>& var,
      const std::string& op_name);

   tbox::Pointer<hier::TimeInterpolateOperator> buildTimeInterpolateOperator(
      const tbox::Pointer<hier::Variable>& var,
      const std::string& op_name);
};

}
}
#endif
