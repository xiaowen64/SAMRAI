/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Common MappedBox operations for MappedBox containers. 
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxContainerUtils
#define included_hier_MappedBoxContainerUtils

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/Connector.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Utilities for performing simple common tasks on a container
 * of MappedBoxes.
 *
 * TODO: Arguments should be re-ordered to the SAMRAI standard, output
 * before input.
 *
 * TODO: There are some very similar methods in this class.  In many cases,
 * one version supports input and output being the same object, but the other
 * does not.  For uniformity, all these methods should support input and
 * output containers being the same object.  It's simple to implement.
 */
class MappedBoxContainerUtils
{

public:


   //@{
   //! @name Changing a plain array of MappedBoxes in place.

   /*!
    * @brief Refine the boxes in a vector of MappedBoxes.
    *
    * @brief[in,out] mapped_box_vector
    *
    * @param[in] ratio Ratio in refinement operation.
    */
   static void
   refineMappedBoxVectorBoxes(
      std::vector<MappedBox>& mapped_box_vector,
      const IntVector& ratio);

   /*!
    * @brief Coarsen the boxes in a vector of MappedBoxes.
    *
    * @brief[in,out] mapped_box_vector
    *
    * @param[in] ratio Ratio in coarsen operation.
    */
   static void
   coarsenMappedBoxVectorBoxes(
      std::vector<MappedBox>& mapped_box_vector,
      const IntVector& ratio);

   /*!
    * @brief Grow the boxes in a vector of MappedBoxes.
    *
    * @param[in,out] mapped_box_vector
    *
    * @param[in] growth Growth amount.
    */
   static void
   growMappedBoxVectorBoxes(
      std::vector<MappedBox>& mapped_box_vector,
      const IntVector& growth);

   //@}



   //@{

   //! @name I/O operations for containers that lack built-in versions.

   /*!
    * @brief Print a vector of MappedBoxes to an output stream.
    *
    * @param[in] mapped_boxes
    *
    * @param[in] output_stream
    *
    * @param[in] left_border
    *
    * @param[in] detail_depth
    */
   static void
   recursivePrintMappedBoxVector(
      const std::vector<MappedBox>& mapped_boxes,
      std::ostream& output_stream = tbox::plog,
      const std::string& border = std::string(),
      int detail_depth = 0);

   //@}



   //@{

   //! @name Some nasty conversions which should eventually go away.

   /**
    * @brief Convert a BoxList to a vector<MappedBox>.
    *
    * Each input box is converted to a MappedBox owned by process 0
    * The LocalIndices are set sequentially starting with 0 for the
    * first input Box.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[out] mapped_box_vector
    *
    * @param[in] block_id
    *
    * @param[in] box_list
    */
   static void
   convertBoxListToMappedBoxVector(
      const BoxList& box_list,
      std::vector<MappedBox>& mapped_box_vector,
      const BlockId& block_id = BlockId::zero());

   /**
    * @brief Convert a vector<MappedBox> to a BoxList.
    *
    * Put the results in the output container.  For flexibility and
    * efficiency, the output container is NOT cleared first, so you
    * may want to clear it before calling this method.
    *
    * @param[in] mapped_box_vector
    * @param[out] box_list
    */
   static void
   convertMappedBoxVectorToBoxList(
      const std::vector<MappedBox>& mapped_box_vector,
      BoxList& box_list);

   //@}

private:

   // Disabled constructor.  No need for objects of this class.
   MappedBoxContainerUtils();


};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxContainerUtils.I"
#endif

#endif  // included_hier_MappedBoxContainerUtils
