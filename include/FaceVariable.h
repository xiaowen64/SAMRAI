//
// File:	FaceVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining face centered variables
//

#ifndef included_pdat_FaceVariable
#define included_pdat_FaceVariable

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
 * Class FaceVariable<DIM> is a templated variable class used to define 
 * face-centered quantities on an AMR mesh.  It is templated on the type 
 * of the underlying data (e.g., double, int, bool, etc.).   Face variable 
 * data is associated with the faces (or sides) of cells.  Face data is 
 * stored in DIM arrays, each of which holds values for faces having the 
 * same normal vector.  For example, a three-dimensional face variable 
 * can be used to create face-centered data arrays over a box 
 * [l0:u0,l1:u1,l2:u2] that can be dimensioned as: 
 * \verbatim

     [ l0 : u0+1 ,
       l1 : u1 ,
       l2 : u2 , d ]   ,

     [ l1 : u1+1 ,
       l2 : u2 ,
       l0 : u0 , d ]   ,

     [ l2 : u2+1 ,
       l0 : u0 ,
       l1 : u1 , d ]   ,

 * \endverbatim
 * for the x, y, and z (or 0, 1, 2) face directions, respectively, and
 * where d is the depth index (i.e., number of values at each face index
 * location).  Note that the array orderings are permuted so that the leading
 * dimension corresponds to the face direction of array.  One- and two-
 * dimensional face variables define storage similarly.  For more information 
 * on indexing and manipulating face patch data objects, see the classes 
 * FaceData<DIM> and FaceGeometry<DIM>.
 *
 * IMPORTANT: The class SideVariable<DIM> and associated classes define
 * the same storage as this face variable class, except that the indices
 * are not permuted in the side data type.
 *
 * @see pdat::FaceData
 * @see pdat::FaceDataFactory
 * @see pdat::FaceGeometry
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class FaceVariable : public hier::Variable<DIM>
{
public:
   /**
    * Create a face variable object having properties specified by the
    * name, depth (i.e., number of data values at each index location),
    * and coarse-fine interface representation.  Default arguments are 
    * provided for the last two.  The default depth is one.  The 
    * fine boundary representation boolean value indicates which values 
    * (either coarse or fine) take precedence during coarsen and refine 
    * operations.  The default state is that fine data values take precedence 
    * on coarse-fine interfaces.  
    */
   FaceVariable(const string &name,
                      int depth = 1,
                      const bool fine_boundary_represents_var = true);

   /**
    * Virtual destructor for face variable objects.
    */
   virtual ~FaceVariable<DIM,TYPE>();

   /**
    * Return a boolean value indicating how data for the face variable will be treated
    * on coarse-fine interfaces.  True (default case set in constructor) indicates
    * that fine patch values take precedence.  False indicates that values on fine patches
    * at a coarse-fine interface should be interpolated from coarser level values.
    */
   bool fineBoundaryRepresentsVariable() const {return d_fine_boundary_represents_var;}

   /**
    * Return true since the face data index space extends beyond the interior of
    * patches.  That is, face data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   bool d_fine_boundary_represents_var;

   FaceVariable(const FaceVariable<DIM,TYPE>&); // not implemented
   void operator=(const FaceVariable<DIM,TYPE>&);	// not implemented
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceVariable.C"
#endif
