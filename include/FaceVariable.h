//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/face/FaceVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Variable class for defining face centered variables
//

#ifndef included_pdat_FaceVariable
#define included_pdat_FaceVariable

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
 * Class FaceVariable<DIM> is a templated variable class used to define
 * face-centered quantities on an AMR mesh.   It is a subclass of
 * hier::Variable and is templated on the type of the underlying data
 * (e.g., double, int, bool, etc.).
 *
 * Note that the indices in the face data arrays are permuted so that
 * the leading index in each array corresponds to the associated face
 * normal coordinate direction. See header file for FaceData<DIM> class 
 * for a more detailed description of the data layout.
 *
 * IMPORTANT: The class SideVariable<DIM> and associated "side data" classes 
 * define the same storage as this face variable class, except that the 
 * individual array indices are not permuted in the side data type.
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
   /*!
    * @brief Create an face-centered variable object with the given name and
    * depth (i.e., number of data values at each edge index location).
    * A default depth of one is provided.   The fine boundary representation
    * boolean argument indicates which values (either coarse or fine) take
    * precedence at coarse-fine mesh boundaries during coarsen and refine
    * operations.  The default is that fine data values take precedence
    * on coarse-fine interfaces.
    */
   FaceVariable(const std::string &name,
                int depth = 1,
                bool fine_boundary_represents_var = true);

   /*!
    * @brief Virtual destructor for face variable objects.
    */
   virtual ~FaceVariable<DIM,TYPE>();

   /*!
    * @brief Return boolean indicating which face data values (coarse
    * or fine) take precedence at coarse-fine mesh interfaces.  The
    * value is set in the constructor.
    */
   bool fineBoundaryRepresentsVariable() const
       {return d_fine_boundary_represents_var;}

   /*!
    * @brief Return true indicating that face data on a patch interior
    * exists on the patch boundary.
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
