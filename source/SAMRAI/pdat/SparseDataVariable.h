/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   SparseDataVariable
 *
 ************************************************************************/
#ifndef included_pdat_SparseDataVariable
#define included_pdat_SparseDataVariable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Variable.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Variable class used to define sparse data on a cell-centered index
 * set.
 *
 * @see SparseData
 * @see SparseDataFactory
 * @see Variable
 */
template <typename BOX_GEOMETRY>
class SparseDataVariable : public hier::Variable
{
public:
   /*!
    * @brief Create a sparse data variable object with the specified name.
    *
    * The creation of the variable creates a "default" SparseDataFactory
    * with ghost width set to zero.
    * 
    * @param [in] dim
    * @param [in] name
    */
   explicit SparseDataVariable(
      const tbox::Dimension& dim,
      const std::string& name,
      const std::vector<std::string>& dbl_attributes,
      const std::vector<std::string>& int_attributes);

   /*!
    * @brief Destructor
    */
   ~SparseDataVariable();

   /*!
    * @brief Returns true.
    *
    * Sparse data quantities will always be treated as cell-centered 
    * quantities for communication purposes.  Note that this is 
    * artificial, since the cell data index space matches the cell-centered
    * index space for AMR patches. Cell data does not live on patch borders,
    * so there is no ambiguity regarding coarse-fine interface values.
    *
    * @see SparseDataVariable::dataLivesOnPatchBorder
    */
   bool
   fineBoundaryRepresentsVariable() const { return true; }

   /*!
    * @brief Returns false.
    * 
    * The sparse data index space matches the cell-centered index space for
    * AMR patches.  Sparse data, therefore, does not live on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const { return false; }

private:
   /*
    * copy c'tor and assignment operator are private to prevent the
    * compiler from generating a default.
    */ 
    SparseDataVariable(
      const SparseDataVariable<BOX_GEOMETRY>& rhs);

    SparseDataVariable<BOX_GEOMETRY>&
    operator=(
      const SparseDataVariable<BOX_GEOMETRY>& rhs);

}; // end class SparseDataVariable.

} // end namespace pdat.
} // end namespace SAMRAI

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SAMRAI/pdat/SparseDataVariable.C"
#endif

#endif
