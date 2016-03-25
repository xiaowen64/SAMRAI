//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/outernode/OuternodeVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Variable class for defining outernode centered variables
//

#ifndef included_pdat_OuternodeVariable
#define included_pdat_OuternodeVariable

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
 * @brief Class OuternodeVariable<DIM> is a templated variable class
 * used to define node-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).
 *
 * Note that the data layout in the outernode data arrays matches the corresponding 
 * array sections provided by the node data implementation.  See header file for 
 * the OuternodeData<DIM> class for a more detailed description of the data layout.
 *
 * @see NodeData<DIM>
 * @see OuternodeData<DIM>
 * @see OuternodeDataFactory<DIM>
 * @see hier::Variable<DIM>
 */

template <int DIM, class TYPE>
class OuternodeVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create an outernode variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.
    *
    * Note that The ghost cell width for all outernode data is currently
    * fixed at zero; this may be changed in the future if needed.
    */
   OuternodeVariable(const std::string &name,
                     int depth = 1);
 
   /*!
    * @brief Virtual destructor for outernode variable objects.
    */
   virtual ~OuternodeVariable<DIM,TYPE>();
 
   /*!
    * @brief Return a boolean true value indicating that fine patch
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}
 
   /*!
    * @brief Return true indicating that outernode data
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   // neither of the following functions are implemented
   OuternodeVariable<DIM,TYPE>(const OuternodeVariable<DIM,TYPE>&);
   void operator=(const OuternodeVariable<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuternodeVariable.C"
#endif


