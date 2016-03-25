//
// File:	OuteredgeVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name$
// Revision:	$Revision: 601 $
// Modified:	$Date: 2005-09-06 11:23:15 -0700 (Tue, 06 Sep 2005) $
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
 * used to define edge-centered quantities on patch boundaries.
 *
 * It is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).  Outeredge variable
 * data is associated with the edges of cells.  However, it differs
 * from the EdgeVariable class in that outeredge quantities reside only
 * on the sides residing on the boundary of a patch.  
 *
 * Outeredge data is stored in DIM*DIM*2 arrays, containing the data for the 
 * patch boundary sides with each of the possible outward pointing normal 
 * directions. Where an outeredge falls on more than one side (patch edges 
 * and corners), the outeredge belongs to the array associated with the 
 * higher dimensional direction. In each of these arrays, memory allocation 
 * is in column-major ordering (e.g., Fortran style) so that the leftmost 
 * index runs fastest in memory.  For example, a three-dimensional outeredge 
 * data object instantiated with a box [l0:u0,l1:u1,l2:u2] allocates 12
 * data (i.e., 3x2 pairs) arrays dimensioned as:
 * \verbatim
 *
 *        (i,j,s) 
 *    X:  (X,Y,[0,1])  [l0:l0,l1:u1,l2+1:u2,d],  [u0:u0,l1:u1,l2+1:u2,d]
 *        (X,Z,[0,1])  [l0:l0,l1+1:u1,l2:u2,d],  [u0:u0,l1+1:u1,l2:u2,d]
 *
 *    Y:  (Y,X,[0,1])  [l0:u0,l1:l1,l2+1:u2,d],  [l0:u0,u1:u1,l2+1:u2,d]
 *        (Y,Z,[0,1])  [l0:u0+1,l1:l1,l2:u2,d],  [l0:u0+1,u1:u1,l2:u2,d]
 *
 *    Z:  (Z,X,[0,1])  [l0:u0,l1:u1+1,l2:l2,d],  [l0:u0,l1:u1+1,u2:u2,d]
 *        (Z,Y,[0,1])  [l0:u0+1,l1:u1,l2:l2,d],  [l0:u0+1,l1:u1,u2:u2,d]
 *
 * \endverbatim
 * where X,Y,and Z can be specified 0, 1, 2, respectively.  
 * One- and two-dimensional edge data arrays are managed similary.  The
 * "j" dimension corresponds with the "axis" of standard EdgeData<DIM>. i.e.
 *
 *    Outeredge box(i,j,s) =  EdgeData<DIM>.getBox(j) 
 *
 * The specific data boxes are constructed based on the value of i and s. 
 * To avoid duplication of outeredge values, the data boxes in lower 
 * dimensions of i are trimmed.  That is, higher dimensions (e.g. Z) take
 * precedent over lower dimensions (e.g. X) if the databox could be defined
 * over both.
 *
 * For more information on indexing and manipulating outeredge patch data 
 * objects, see the classes OuteredgeData<DIM> and OuteredgeGeometry<DIM>.
 *
 * @see pdat::EdgeData
 * @see pdat::OuteredgeData
 * @see pdat::OuteredgeDataFactory
 * @see hier::Variable<DIM>
 */

template <int DIM, class TYPE>
class OuteredgeVariable : public hier::Variable<DIM>
{
public:
   /*!
    * @brief Create an outeredge variable object having properties
    * specified by the name and depth (i.e., number of data values
    * at each index location).  The default depth is one.   The ghost 
    * cell width for all outeredge data is currently fixed at zero; 
    * this may be changed in the future if needed.
    */
   OuteredgeVariable(const string &name, 
                     int depth = 1);

   /*!
    * @brief Virtual destructor for outeredge variable objects.
    */
   virtual ~OuteredgeVariable<DIM,TYPE>();

   /*!
    * @brief Return a boolean true value indicating that fine patch values take precedence
    * on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /*!
    * Return true since the edge data index space (and hence the outeredge data index
    * space) extends beyond the interior of patches.  That is, outeredge data lives
    * on patch borders.
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


