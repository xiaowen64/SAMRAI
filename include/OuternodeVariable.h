//
// File:	OuternodeVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Variable class for defining outernode centered variables
//

#ifndef included_pdat_OuternodeVariable
#define included_pdat_OuternodeVariable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class OuternodeVariable<DIM> is a templated variable class
 * used to define node-centered quantities on patch boundaries.
 *
 * It is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).  Outernode variable
 * data is associated with the nodes of cells.  However, it differs
 * from the NodeVariable<DIM> class in that outernode quantities reside only
 * on the sides residing on the boundary of a patch.  Outernode data is
 * stored in 2*DIM arrays, each of which holds values for sides having the
 * same outward normal vector.
 * Where an outernode falls on more than one side (patch edges and corners),
 * the outernode belongs to the higher dimensional direction.
 * For example, a three-dimensional outernode
 * variable can be used to create node-centered data arrays over a box
 * [l0:u0,l1:u1,l2:u2] that can be dimensioned as:
 * \verbatim
 *
 *    [ l1+1 : u1-1 ,
 *      l2+1 : u2-1 , d ]   ,
 *
 *    [ l0   : u0   ,
 *      l2+1 : u2-1 , d ]   ,
 *
 *    [ l0   : u0   ,
 *      l1   : u1   , d ]   ,
 *
 * \endverbatim
 * for the upper and lower x, y, and z (or 0, 1, 2) face directions,
 * respectively, and where d is the depth index (i.e., number of values at
 * each side index location).  One- and two-dimensional 
 * outernode variables define storage similarly.  For more information on 
 * indexing and manipulating outernode patch data objects, see the classes 
 * OuternodeData<DIM> and OuternodeGeometry<DIM>.
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
    * at each index location).  The default depth is one.   The ghost 
    * cell width for all outernode data is currently fixed at zero; 
    * this may be changed in the future if needed.
    */
   OuternodeVariable(const string &name, 
                     int depth = 1);

   /*!
    * @brief Virtual destructor for outernode variable objects.
    */
   virtual ~OuternodeVariable<DIM,TYPE>();

   /*!
    * @brief Return a boolean true value indicating that fine patch values 
    * take precedence on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /*!
    * Return true since the node data index space (and hence the outernode 
    * data index space) extends beyond the interior of patches.  That is, 
    * outernode data lives on patch borders.
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


