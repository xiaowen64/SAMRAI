//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/cell/CellVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Variable class for defining cell centered variables
//

#ifndef included_pdat_CellVariable
#define included_pdat_CellVariable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_String
#include <string>
#define included_String
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif

namespace SAMRAI {
    namespace pdat {

/*!
 * Class CellVariable<DIM> is a templated variable class used to define 
 * cell-centered quantities on an AMR mesh.   It is a subclass of 
 * hier::Variable and is templated on the type of the underlying data 
 * (e.g., double, int, bool, etc.).
 *
 * See header file for CellData<DIM> class for a more detailed
 * description of the data layout. 
 *
 * @see pdat::CellData
 * @see pdat::CellDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class CellVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create a cell-centered variable object with the given name and
    * depth (i.e., number of data values at each cell index location).
    * A default depth of one is provided. 
    */
   CellVariable(const std::string &name,
                int depth = 1);

   /*!
    * @brief Virtual destructor for cell variable objects.
    */
   virtual ~CellVariable<DIM,TYPE>();

   /*!
    * @brief Return true indicating that cell data quantities will always 
    * be treated as though fine values take precedence on coarse-fine 
    * interfaces.  Note that this is really artificial since the cell 
    * data index space matches the cell-centered index space for AMR 
    * patches.  However, some value must be supplied for communication 
    * operations.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /*!
    * @brief Return false indicating that cell data on a patch interior
    * does not exist on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return false;}

private:
   CellVariable(const CellVariable<DIM,TYPE>&);// not implemented
   void operator=(const CellVariable<DIM,TYPE>&);	// not implemented
};


}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellVariable.C"
#endif
