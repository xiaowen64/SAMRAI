//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outerface/OuterfaceVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Variable class for defining outerface centered variables
//

#ifndef included_pdat_OuterfaceVariable
#define included_pdat_OuterfaceVariable

#include "SAMRAI_config.h"
#include "tbox/Complex.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "Variable.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class OuterfaceVariable<DIM> is a templated variable class
 * used to define face-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).
 *
 * Note that the data layout in the outerface data arrays matches the corresponding 
 * array sections provided by the face data implementation. See header file for 
 * the OuterfaceData<DIM> class for a more detailed description of the data layout.
 *
 * IMPORTANT: The class OutersideVariable<DIM> and associated "outerside 
 * data" classes define the same storage as this outerface variable class, 
 * except that the individual array indices are not permuted in the outerside 
 * data type.
 *
 * @see pdat::FaceData
 * @see pdat::OuterfaceData
 * @see pdat::OuterfaceDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class OuterfaceVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create an outerface variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.
    *
    * Note that The ghost cell width for all outerface data is currently
    * fixed at zero; this may be changed in the future if needed.
    */
   OuterfaceVariable(const std::string &name,
                     int depth = 1);

   /*!
    * @brief Virtual destructor for outerface variable objects.
    */
   virtual ~OuterfaceVariable<DIM,TYPE>();
 
   /*!
    * @brief Return a boolean true value indicating that fine patch
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}
 
   /*!
    * @brief Return true indicating that outerface data
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   // neither of the following functions are implemented
   OuterfaceVariable(const OuterfaceVariable<DIM,TYPE>&);
   void operator=(const OuterfaceVariable<DIM,TYPE>&);
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuterfaceVariable.C"
#endif
