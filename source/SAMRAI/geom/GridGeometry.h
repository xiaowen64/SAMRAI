/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management in AMR hierarchy
 *
 ************************************************************************/

#ifndef included_geom_GridGeometry
#define included_geom_GridGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include "boost/shared_ptr.hpp"

namespace SAMRAI {
namespace geom {

/*!
 * @brief Class GridGeometry manages the index space that determines the
 * extent of the coarse-level domain of a SAMRAI hierarchy.
 *
 * A GridGeometry object can be directly constructed with a consistent
 * state, or it may be used as a base class to derive child classes that manage
 * particular grid types (%e.g., Cartesian, cylindrical, etc.).  Direct
 * construction of a GridGeometry object is recommended for multiblock
 * problems and other problems where the physical locations of mesh
 * coordinates are managed by user code.
 *
 * @see hier::BaseGridGeometry
 */

class GridGeometry:
   public hier::BaseGridGeometry
{
   friend class TransferOperatorRegistry;

public:
   /*!
    * @brief Construct a new GridGeometry object and initialize from
    * input.
    *
    * This constructor for GridGeometry initializes data members
    * based on parameters read from the specified input database.
    *
    * This constructor is intended for use when directly constructing a
    * GridGeometry without using a derived child class.  The object will
    * contain all index space grid information for a mesh, but nothing about
    * the physical coordinates of the mesh.
    *
    * @note
    * @b Errors: passing in a null database pointer or an empty string
    * will result in an unrecoverable assertion.
    *
    * @param[in]  dim
    * @param[in]  object_name
    * @param[in]  input_db
    */
   GridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const boost::shared_ptr<tbox::Database>& input_db);

   /*!
    * @brief Construct a new GridGeometry object based on arguments.
    *
    * This constructor creates a new GridGeometry object based on the
    * arguments, rather than relying on input or restart data.
    *
    * @param[in]  object_name
    * @param[in]  domain      Each element of the array describes the index
    *                         space for a block.
    */
   GridGeometry(
      const std::string& object_name,
      const hier::BoxContainer& domain);

   /*!
    * @brief Construct a new coarsened/refined GridGeometry object with the
    * supplied domain.
    *
    * This method is intended to be called only by boost::make_shared from the
    * make[Coarsened, Refined]GridGeometry methods to make a coarsened or
    * refined version of a given GridGeometry.
    *
    * @param[in]  object_name The same name as the uncoarsened/unrefined grid
    *                         geometry.
    * @param[in]  domain The coarsened/refined domain.
    * @param[in]  op_reg The same operator registry as the
    *                    uncoarsened/unrefined grid geometry.
    */
   GridGeometry(
      const std::string& object_name,
      const hier::BoxContainer& domain,
      const boost::shared_ptr<hier::TransferOperatorRegistry>& op_reg);

   /*!
    * @brief Virtual destructor
    */
   virtual ~GridGeometry();

   /*!
    * @brief Create a pointer to a refined version of this grid geometry
    *        object.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in]     fine_geom_name std::string name of the geometry object
    * @param[in]     refine_ratio the refinement ratio.
    *
    * @return The pointer to the grid geometry object.
    */
   virtual boost::shared_ptr<hier::BaseGridGeometry>
   makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector& refine_ratio) const;

   /*!
    * @brief Create a pointer to a coarsened version of this grid geometry
    *        object.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in]     coarse_geom_name std::string name of the geometry object
    * @param[in]     coarsen_ratio the coasening ratio
    *
    * @return The pointer to a coarsened version of this grid geometry object.
    */
   virtual boost::shared_ptr<hier::BaseGridGeometry>
   makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector& coarsen_ratio) const;

protected:
   /*!
    * @brief Build operators appropriate for a GridGeometry.
    */
   virtual void
   buildOperators();
};

}
}

#endif
