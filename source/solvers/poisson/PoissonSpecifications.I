/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/PoissonSpecifications.I $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Specifications for the scalar Poisson equation
 */

#include "tbox/Utilities.h"
#include "tbox/Pointer.h"


#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
    namespace solv {

/*
*******************************************************************
* Default constructor                                             *
*******************************************************************
*/

inline
PoissonSpecifications::PoissonSpecifications(
   const std::string &object_name
   ) : d_object_name(object_name),
       d_D_id(-1),
       d_D_constant(1.0),
       d_C_zero(true),
       d_C_id(-1),
       d_C_constant(0.0) {
   return;
}


/*
*******************************************************************
* Copy constructor                                                *
*******************************************************************
*/

inline
PoissonSpecifications::PoissonSpecifications(
   const std::string &object_name,
   const PoissonSpecifications &r
   ) : d_object_name(object_name),
       d_D_id(r.d_D_id),
       d_D_constant(r.d_D_constant),
       d_C_zero(r.d_C_zero),
       d_C_id(r.d_C_id),
       d_C_constant(r.d_C_constant) {
   return;
}

/*
*******************************************************************
* Destructor (does nothing).                                      *
*******************************************************************
*/
inline
PoissonSpecifications::~PoissonSpecifications() {
   return;
}


/*
*******************************************************************
* Assignment operator                                             *
*******************************************************************
*/

inline
const PoissonSpecifications
&PoissonSpecifications::operator=(
   const PoissonSpecifications &r ) {
   d_D_id = r.d_D_id;
   d_D_constant = r.d_D_constant;
   d_C_zero = r.d_C_zero;
   d_C_id = r.d_C_id;
   d_C_constant = r.d_C_constant;
   return *this;
}

/*
*******************************************************************
* Set the patch data index for variable D.                        *
*******************************************************************
*/
inline
void PoissonSpecifications::setDPatchDataId( int id ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( id < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_D_id = id;
   d_D_constant = 0.0;
   return;
}

/*
*******************************************************************
* Set the constant value variable D.                              *
*******************************************************************
*/
inline
void PoissonSpecifications::setDConstant( double constant ) {
   d_D_id = -1;
   d_D_constant = constant;
   return;
}

/*
*******************************************************************
* Whether D is variable.                                          *
*******************************************************************
*/
inline
bool PoissonSpecifications::dIsVariable() const {
   return d_D_id != -1;
}

/*
*******************************************************************
* Whether D is constant.                                          *
*******************************************************************
*/
inline
bool PoissonSpecifications::dIsConstant() const {
   return d_D_id == -1;
}

/*
*******************************************************************
* Get the patch data index for variable D.                        *
*******************************************************************
*/
inline
int PoissonSpecifications::getDPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_D_id == -1 ) {
      TBOX_ERROR(d_object_name << ": D not prepresented by a patch data.\n");
   }
#endif
   return d_D_id;
}

/*
*******************************************************************
* Get the constant D value.                                       *
*******************************************************************
*/
inline
double PoissonSpecifications::getDConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_D_id != -1 ) {
      TBOX_ERROR(d_object_name << ": D not prepresented by a constant.\n");
   }
#endif
   return d_D_constant;
}

/*
*******************************************************************
* Set the constant value variable C.                              *
*******************************************************************
*/
inline
void PoissonSpecifications::setCPatchDataId( int id ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( id < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_C_zero = false;
   d_C_id = id;
   d_C_constant = 0.0;
   return;
}

/*
*******************************************************************
* Set the patch data index for variable C.                        *
*******************************************************************
*/
inline
void PoissonSpecifications::setCConstant( double constant ) {
   d_C_zero = false;
   d_C_id = -1;
   d_C_constant = constant;
   return;
}

/*
*******************************************************************
* Set the value of C to zero.                                     *
*******************************************************************
*/
inline
void PoissonSpecifications::setCZero() {
   d_C_zero = true;
   d_C_id = -1;
   d_C_constant = 0.0;
   return;
}

/*
*******************************************************************
* Whether C is variable.                                          *
*******************************************************************
*/
inline
bool PoissonSpecifications::cIsVariable() const {
   return d_C_id != -1;
}

/*
*******************************************************************
* Whether C is zero.                                              *
*******************************************************************
*/
inline
bool PoissonSpecifications::cIsZero() const {
   return d_C_zero;
}

/*
*******************************************************************
* Whether C is constant.                                          *
*******************************************************************
*/
inline
bool PoissonSpecifications::cIsConstant() const {
   return !d_C_zero && (d_C_id == -1);
}

/*
*******************************************************************
* Get the patch data index for variable C.                        *
*******************************************************************
*/
inline
int PoissonSpecifications::getCPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_C_id == -1 ) {
      TBOX_ERROR(d_object_name << ": C not prepresented by a an index.\n");
   }
#endif
   return d_C_id;
}

/*
*******************************************************************
* Get the constant C value.                                       *
*******************************************************************
*/
inline
double PoissonSpecifications::getCConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_C_id != -1 || d_C_zero ) {
     TBOX_ERROR(d_object_name << ": C is not prepresented by a constant.\n");
   }
#endif
   return d_C_constant;
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif

} // namespace SAMRAI
}
