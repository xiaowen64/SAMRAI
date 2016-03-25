//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outeredge/OuteredgeVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Variable class for defining outeredge centered variables
//

#ifndef included_pdat_OuteredgeVariable
#define included_pdat_OuteredgeVariable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class OuteredgeVariable<DIM> is a templated variable class
 * used to define edge-centered data quantities only on patch boundaries.
 * It is a subclass of hier::Variable and is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).  
 *
 * Note that the data layout in the outeredge data arrays matches the corresponding 
 * array sections provided by the edge data implementation.  See header file for 
 * the OuteredgeData<DIM> class for a more detailed description of the data layout. 
 *
 * @see pdat::EdgeData
 * @see pdat::OuteredgeData
 * @see pdat::OuteredgeDataFactory
 * @see hier::Variable
 */

template <int DIM, class TYPE>
class OuteredgeVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create an outeredge variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.   
    * 
    * Note that The ghost cell width for all outeredge data is currently 
    * fixed at zero; this may be changed in the future if needed.
    */
   OuteredgeVariable(const std::string &name, 
                     int depth = 1);

   /*!
    * @brief Virtual destructor for outeredge variable objects.
    */
   virtual ~OuteredgeVariable<DIM,TYPE>();

   /*!
    * @brief Return a boolean true value indicating that fine patch 
    * values take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /*!
    * @brief Return true indicating that outeredge data 
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   // neither of the following functions are implemented
   OuteredgeVariable(const OuteredgeVariable<DIM,TYPE>&);
   void operator=(const OuteredgeVariable<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuteredgeVariable.C"
#endif


