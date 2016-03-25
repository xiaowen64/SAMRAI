//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/node/NodeVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	hier::Variable class for defining node centered variables
//

#ifndef included_pdat_NodeVariable
#define included_pdat_NodeVariable

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
 * Class NodeVariable<DIM> is a templated variable class used to define
 * node-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * See header file for NodeData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see pdat::NodeData
 * @see pdat::NodeDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class NodeVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create a cell-centered variable object with the given name and
    * depth (i.e., number of data values at each cell index location).
    * A default depth of one is provided.  The fine boundary representation 
    * boolean argument indicates which values (either coarse or fine) take 
    * precedence at coarse-fine mesh boundaries during coarsen and refine 
    * operations.  The default is that fine data values take precedence 
    * on coarse-fine interfaces.
    */
   NodeVariable(const std::string &name,
                int depth = 1, 
                bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for node variable objects.
    */
   virtual ~NodeVariable<DIM,TYPE>();

   /*!
    * @brief Return boolean indicating which node data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const 
       {return d_fine_boundary_represents_var;}

   /*!
    * @brief Return true indicating that node data on a patch interior
    * exists on the patch boundary.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   bool d_fine_boundary_represents_var;

   NodeVariable(const NodeVariable<DIM,TYPE>&);// not implemented
   void operator=(const NodeVariable<DIM,TYPE>&);	// not implemented
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NodeVariable.C"
#endif
