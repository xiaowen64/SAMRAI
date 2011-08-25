/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Pointer.h"

namespace SAMRAI {
namespace solv {

/*
 *******************************************************************
 * Default constructor                                             *
 *******************************************************************
 */

SAMRAI_INLINE_KEYWORD
PoissonSpecifications::PoissonSpecifications(
   const std::string& object_name):d_object_name(object_name),
   d_D_id(-1),
   d_D_constant(1.0),
   d_C_zero(true),
   d_C_id(-1),
   d_C_constant(0.0) {
}

/*
 *******************************************************************
 * Copy constructor                                                *
 *******************************************************************
 */

SAMRAI_INLINE_KEYWORD
PoissonSpecifications::PoissonSpecifications(
   const std::string& object_name,
   const PoissonSpecifications& r):d_object_name(object_name),
   d_D_id(r.d_D_id),
   d_D_constant(r.d_D_constant),
   d_C_zero(r.d_C_zero),
   d_C_id(r.d_C_id),
   d_C_constant(r.d_C_constant) {
}

/*
 *******************************************************************
 * Destructor (does nothing).                                      *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
PoissonSpecifications::~PoissonSpecifications() {
}

/*
 *******************************************************************
 * Assignment operator                                             *
 *******************************************************************
 */

SAMRAI_INLINE_KEYWORD
const PoissonSpecifications
& PoissonSpecifications::operator = (
   const PoissonSpecifications& r) {
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
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setDPatchDataId(
   int id) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (id < 0) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_D_id = id;
   d_D_constant = 0.0;
}

/*
 *******************************************************************
 * Set the constant value variable D.                              *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setDConstant(
   double constant) {
   d_D_id = -1;
   d_D_constant = constant;
}

/*
 *******************************************************************
 * Whether D is variable.                                          *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::dIsVariable() const {
   return d_D_id != -1;
}

/*
 *******************************************************************
 * Whether D is constant.                                          *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::dIsConstant() const {
   return d_D_id == -1;
}

/*
 *******************************************************************
 * Get the patch data index for variable D.                        *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
int PoissonSpecifications::getDPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_D_id == -1) {
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
SAMRAI_INLINE_KEYWORD
double PoissonSpecifications::getDConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_D_id != -1) {
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
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setCPatchDataId(
   int id) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (id < 0) {
      TBOX_ERROR(d_object_name << ": Invalid patch data id.\n");
   }
#endif
   d_C_zero = false;
   d_C_id = id;
   d_C_constant = 0.0;
}

/*
 *******************************************************************
 * Set the patch data index for variable C.                        *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setCConstant(
   double constant) {
   d_C_zero = false;
   d_C_id = -1;
   d_C_constant = constant;
}

/*
 *******************************************************************
 * Set the value of C to zero.                                     *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
void PoissonSpecifications::setCZero() {
   d_C_zero = true;
   d_C_id = -1;
   d_C_constant = 0.0;
}

/*
 *******************************************************************
 * Whether C is variable.                                          *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::cIsVariable() const {
   return d_C_id != -1;
}

/*
 *******************************************************************
 * Whether C is zero.                                              *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::cIsZero() const {
   return d_C_zero;
}

/*
 *******************************************************************
 * Whether C is constant.                                          *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool PoissonSpecifications::cIsConstant() const {
   return !d_C_zero && (d_C_id == -1);
}

/*
 *******************************************************************
 * Get the patch data index for variable C.                        *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
int PoissonSpecifications::getCPatchDataId() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_C_id == -1) {
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
SAMRAI_INLINE_KEYWORD
double PoissonSpecifications::getCConstant() const {
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_C_id != -1 || d_C_zero) {
      TBOX_ERROR(d_object_name << ": C is not prepresented by a constant.\n");
   }
#endif
   return d_C_constant;
}

/*
 *******************************************************************
 * Get the object's name.                                          *
 *******************************************************************
 */
SAMRAI_INLINE_KEYWORD
const std::string& PoissonSpecifications::getObjectName() const {
   return d_object_name;
}

} // namespace SAMRAI
}
