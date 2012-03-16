/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevelStatistics
#define included_hier_BoxLevelStatistics

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxLevel.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief A utility for writing out various statistics of MappedBoxes.
 */
class BoxLevelStatistics
{

public:

   /*!
    * @brief Constructor.
    *
    * @param[i] box_level
    */
   explicit BoxLevelStatistics(
      const BoxLevel &box_level);

   /*!
    * @brief Print out local and globally reduced statistics on the
    * Boxes.
    *
    * @param[in,out] os The output stream
    *
    * @param[in] border A string to print at the start of every line
    * in the output.
    */
   void
   printBoxStats(
      std::ostream& os,
      const std::string& border) const;


private:

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief Indices for statistical quantites.
    */
   enum { HAS_ANY_BOX,
          NUMBER_OF_CELLS,
          NUMBER_OF_BOXES,
          LARGEST_BOX_SIZE,
          SMALLEST_BOX_SIZE,
          LONGEST_BOX_LEN,
          SHORTEST_BOX_LEN,
          LARGEST_ASPECT_RATIO,
          SMALLEST_ASPECT_RATIO,
          SUM_SURFACE_AREA,
          SUM_NORMALIZED_SURFACE_AREA,
          NUMBER_OF_QUANTITIES };

   /*
    * @brief StatisticalQuantities to compute the min, avg and max for.
    *
    * These quantities will be computed locally on each process and
    * globally reduced.  Not all of these quantities are floating
    * points but all are represented as such.
    */
   struct StatisticalQuantities {
      StatisticalQuantities();
      double d_values[NUMBER_OF_QUANTITIES];
   };

   void
   computeLocalBoxLevelStatistics(
      StatisticalQuantities &sq) const;

   //! @brief BoxLevel to compute statistics for.
   const BoxLevel &d_box_level;

   /*!
    * @brief Names of the quantities in StatisticalQuantities.
    */
   static std::string s_quantity_names[NUMBER_OF_QUANTITIES];

   /*!
    * @brief Longest length in s_quantity_names.
    */
   static int s_longest_length;

   static tbox::StartupShutdownManager::Handler s_initialize_finalize_handler;


};

}
}

#endif  // included_hier_BoxLevelStatistics
