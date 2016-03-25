//
// File:	OuterfaceVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining outerface centered variables
//

#ifndef included_pdat_OuterfaceVariable
#define included_pdat_OuterfaceVariable

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

/**
 * Class OuterfaceVariable<DIM> is a templated variable class used to define
 * face-centered quantities on an AMR mesh.  It is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).  Outerface variable
 * data is associated with the faces (or sides) of cells.  However, it differs
 * from the FaceVariable<DIM> class in that outerface quantities reside only
 * on the faces residing on the boundary of a patch.  Outerface data is
 * stored in DIM arrays, each of which holds values for faces having the
 * same outward normal vector.  For example, a three-dimensional outerface
 * variable can be used to create face-centered data arrays over a box
 * [l0:u0,l1:u1,l2:u2] that can be dimensioned as:
 * \verbatim

     [ l1 : u1 ,
       l2 : u2 , d ]   ,

     [ l2 : u2 ,
       l0 : u0 , d ]   ,

     [ l0 : u0 ,
       l1 : u1 , d ]   ,

 * \endverbatim
 * for the upper and lower x, y, and z (or 0, 1, 2) face directions, 
 * respectively, and where d is the depth index (i.e., number of values at 
 * each face index location).  Note that the array orderings are permuted 
 * to match the conventions of the FaceVariable<DIM> class.  One- and two-
 * dimensional outerface variables define storage similarly.  For more 
 * information on indexing and manipulating outerface patch data objects, 
 * see the classes OuterfaceData<DIM> and OuterfaceGeometry<DIM>.
 *
 * IMPORTANT: The class OutersideVariable<DIM> and associated classes define
 * the same storage as this outerface variable class, except that the indices
 * are not permuted in the outerside data type.  Also, outerface and outerside 
 * data classes are intended to interact with their face-centered and 
 * side-centered data counterparts, respectively.  Mixing types, while 
 * allowed, is discouraged to prevent undesirable behavior.
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
   /**
    * Create an outerface variable object having properties specified by the
    * name and depth (i.e., number of data values at each index location).
    * The default depth is one.   The ghost cell width for all outerface
    * data is currently fixed at zero; this may be changed in the future 
    * if needed.
    */
   OuterfaceVariable(const string &name,
                           int depth = 1);

   /**
    * Virtual destructor for outerface variable objects.
    */
   virtual ~OuterfaceVariable<DIM,TYPE>();

   /**
    * Return a boolean true value indicating that fine patch values take precedence 
    * on coarse-fine interfaces. 
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return true since the face data index space (and hence the outerface data index
    * space) extends beyond the interior of patches.  That is, outerface data lives 
    * on patch borders.
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
