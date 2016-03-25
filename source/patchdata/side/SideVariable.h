//
// File:	SideVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining side centered variables
//

#ifndef included_pdat_SideVariable
#define included_pdat_SideVariable

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
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class SideVariable<DIM> is a templated variable class used to define
 * side-centered quantities on an AMR mesh.  It is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).   Side variable
 * data is associated with the sides (or faces) of cells.  Side data is
 * stored in DIM arrays, each of which holds values for sides having the
 * same normal vector.  For example, a three-dimensional side variable
 * can be used to create side-centered data arrays over a box
 * [l0:u0,l1:u1,l2:u2] that can be dimensioned as:
 * \verbatim

     [ l0 : u0+1 ,
       l1 : u1 ,
       l2 : u2 , d ]   ,

     [ l0 : u0 ,
       l1 : u1+1 ,
       l2 : u2 , d ]   ,

     [ l0 : u0 ,
       l1 : u1 ,
       l2 : u2+1 , d ]   ,

 * \endverbatim
 * for the x, y, and z (or 0, 1, 2) side directions, respectively, and
 * where d is the depth index (i.e., number of values at each side index
 * location).  One- and two- dimensional side variables define storage 
 * similarly.  For more information on indexing and manipulating side 
 * patch data objects, see the classes SideData<DIM> and SideGeometry<DIM>. 
 *
 * IMPORTANT: The class FaceVariable<DIM> and associated classes define
 * the same storage as this side variable class, except that the indices
 * are permuted in the face data type.
 * 
 * Note that it is possible to create a side variable to allocate and 
 * manage data for cell sides associated with a single coordinate direction
 * only.  See the constructor for more information.
 *
 * @see pdat::SideData
 * @see pdat::SideDataFactory
 * @see pdat::SideGeometry
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class SideVariable : public hier::Variable<DIM>
{
public:
   /**
    * Create a side variable object having properties specified by the
    * name, depth (i.e., number of data values at each index location), 
    * coarse-fine interface representation, and coordinate direction 
    * information.  Default arguments are provided for the last three.
    * The default depth is one.  The fine boundary representation boolean
    * indicates which values (either coarse or fine) take precedence 
    * during coarsen and refine operations.  The default state is that 
    * fine data values take precedence on coarse-fine interfaces.  
    * The default data allocation scheme is that side data will 
    * be allocated for all coordinate directions (i.e., -1).  If this is
    * desired, then the direction argument may be omitted.   If an integer
    * direction argument is specified, the only data for that direction
    * will be maintained and managed for this variable (if not -1). 
    */
   SideVariable(const string &name,
                      int depth = 1,
                      bool fine_boundary_represents_var = true,
                      int direction = -1);

   /**
    * Virtual destructor for side variable objects.
    */
   virtual ~SideVariable<DIM,TYPE>();

   /**
    * Return constant reference to vector describing which coordinate
    * directions have data associated with this side variable object.
    * A vector entry of zero indicates that there is no data array
    * allocated for the corresponding coordinate direction for side data
    * created via this side variable object.  A non-zero value indicates
    * that a valid data array will be allocated for that coordinate
    * direction.
    */
   const hier::IntVector<DIM>& getDirectionVector() const;

   /**
    * Return a boolean value indicating how data for the side variable will be treated
    * on coarse-fine interfaces.  True (default case set in constructor) indicates
    * that fine patch values take precedence.  False indicates that values on fine patches
    * at a coarse-fine interface should be interpolated from coarser level values.
    */
   bool fineBoundaryRepresentsVariable() const {return d_fine_boundary_represents_var;}

   /**
    * Return true since the side data index space extends beyond the interior of
    * patches.  That is, side data lives on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   bool d_fine_boundary_represents_var;
   hier::IntVector<DIM> d_directions;
   
   SideVariable(const SideVariable<DIM,TYPE>&);  // not implemented
   void operator=(const SideVariable<DIM,TYPE>&);	// not implemented

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SideVariable.C"
#endif
