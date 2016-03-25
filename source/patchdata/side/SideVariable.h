//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/side/SideVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Variable class for defining side centered variables
//

#ifndef included_pdat_SideVariable
#define included_pdat_SideVariable

#include "SAMRAI_config.h"
#include "tbox/Complex.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "IntVector.h"
#include "Variable.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * Class SideVariable<DIM> is a templated variable class used to define
 * side-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).  See header file for SideData<DIM> class
 * for a more detailed description of the data layout.
 *
 * Note that it is possible to create a side data object for managing
 * data at cell sides associated with a single coordinate direction only.
 * See the constructor for more information. 
 *
 * IMPORTANT: The class FaceVariable<DIM> and associated "face data" classes
 * define the same storage as this side variable class, except that the
 * individual array indices are permuted in the face data type.
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
   /*!
    * @brief Create an side-centered variable object with the given name and
    * depth (i.e., number of data values at each edge index location).
    * A default depth of one is provided.   The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine
    * operations.  The default is that fine data values take precedence
    * on coarse-fine interfaces.  
    *
    * The default data allocation scheme is that side data will
    * be allocated for all side normal coordinate directions.  If this is
    * desired, then the direction argument may be omitted.   If an integer
    * direction argument is specified, the only data for the given 
    * side normal direction will be maintained and managed for this variable.
    */
   SideVariable(const std::string &name,
                int depth = 1,
                bool fine_boundary_represents_var = true,
                int direction = -1);

   /*!
    * @brief Virtual destructor for side variable objects.
    */
   virtual ~SideVariable<DIM,TYPE>();

   /*!
    * @brief Return constant reference to vector describing which coordinate
    * directions have data associated with this side variable object.
    *
    * A vector entry of zero indicates that there is no data array
    * allocated for the corresponding coordinate direction for side data
    * created via this side variable object.  A non-zero value indicates
    * that a valid data array will be allocated for that coordinate
    * direction.
    */
   const hier::IntVector<DIM>& getDirectionVector() const;

   /*!
    * @brief Return boolean indicating which side data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
       {return d_fine_boundary_represents_var;}

   /*!
    * @brief Return true indicating that side data on a patch interior
    * exists on the patch boundary.
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
