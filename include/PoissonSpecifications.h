/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/solvers/poisson/PoissonSpecifications.h $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1845 $
 * Modified:    $LastChangedDate: 2008-01-10 14:47:51 -0800 (Thu, 10 Jan 2008) $
 * Description: Specifications for the scalar Poisson equation
 */

#ifndef included_solv_PoissonSpecifications
#define included_solv_PoissonSpecifications

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_SAMRAI_String
#include <string>
#define included_String
#endif


namespace SAMRAI {
    namespace solv {

/*!
 * @brief Light class holding specifications for cell-centered
 * implementation of the scalar Poisson equation.
 *
 * The scalar Poisson equation is
 * @f$ \nabla ( D \nabla u ) + C u = f @f$,
 * where C is a scalar field, D is the diffusion coefficient.
 * and u and f are scalar quantities.
 *
 * This class describes the things you can set: C, D.
 *
 * Note that the storage and alignment of u, f, C and D depend
 * on the implementation of the solver.  For example, if the
 * solver is cell centered, u, f and C are cell-centered while
 * D is side-centered.
 */

class PoissonSpecifications
{
public:

   /*!
    * @brief Constructor.
    *
    * Sets the specifications to their default state:
    * - C is zero
    * - D is uniformly 1
    *
    * @param object_name Name of object.
    */
   PoissonSpecifications( const std::string &object_name );

   /*!
    * @brief Copy constructor.
    */
   PoissonSpecifications(
      const std::string &object_name,
      const PoissonSpecifications &r );

   /*!
    * @brief Destructor (does nothing).
    */
   ~PoissonSpecifications();

   /*!
    * @brief Assignment operator
    *
    * Assign everything except name.
    */
   const PoissonSpecifications &operator=(
      const PoissonSpecifications &r );

   /*!
    * @brief Print out class data.
    */
   virtual void printClassData( std::ostream &stream ) const;

   //@{
   //! @name Functions for setting and getting D

   /*!
    * @brief Set the patch data index for variable D.
    *
    * In addition, disregard any previous value
    * specified by setDConstant().
    */
   void setDPatchDataId( int id );

   /*!
    * @brief Set the constant value variable D.
    *
    * In addition, disregard any previous patch data index
    * specified by setDPatchDataId().
    */
   void setDConstant( double constant );

   /*!
    * @brief Whether D is variable (described by a patch data id).
    *
    * @return True if D is variable, described by the patch data
    *         id given in setCPatchDataId().
    */
   bool dIsVariable() const;

   /*!
    * @brief Whether D is constant.
    *
    * @return True if D is constant, as specified by setCConstant().
    */
   bool dIsConstant() const;

   /*!
    * @brief Get D's patch data id
    *
    * Error if D is not represented by a patch data id.
    *
    * @return D's id
    */
   int getDPatchDataId() const;

   /*!
    * @brief Get D constant value
    *
    * Error if D is not represented by a constant.
    *
    * @return D's constant value
    */
   double getDConstant() const;

   //@}

   //@{
   //! @name Functions for setting and getting C

   /*!
    * @brief Set the patch data index for C.
    *
    * In addition, disregard any previous values
    * specified by setCConstant() or setCZero().
    */
   void setCPatchDataId( int id );

   /*!
    * @brief Set C to a constant.
    *
    * In addition, disregard any previous value
    * specified by setCPatchDataId() or setCZero().
    *
    * If you want to set C to zero, use setCZero() instead.
    * This allows solvers to take advantage of fact C is absent.
    */
   void setCConstant( double constant );

   /*!
    * @brief Set the value of C to zero.
    *
    * In addition, disregard any previous patch data index
    * specified by setCPatchDataId() and any previous constant
    * specified by setCConstant().
    */
   void setCZero();

   /*!
    * @brief Whether C is variable (described by a patch data id).
    *
    * @return True if C is variable, described by the patch data
    *         id given in setCPatchDataId().
    */
   bool cIsVariable() const;

   /*!
    * @brief Whether C is zero.
    *
    * As it pertains to what this function returns,
    * C is zero @em only by calling setCZero().
    * Calling setCConstant() does @em not make C zero,
    * even if you pass in the value of zero.
    *
    * @return True if C is exactly zero, as set by setCZero().
    */
   bool cIsZero() const;

   /*!
    * @brief Whether C is constant.
    *
    * As it pertains to what this function returns,
    * C is constant @em only by calling setCConstant().
    * Calling setCZero() does @em not make C a constant.
    *
    * @return True if C is constant, as specified by setCConstant().
    */
   bool cIsConstant() const;

   /*!
    * @brief Get C's patch data id
    *
    * Error if C is not represented by a patch data id.
    *
    * @return C's patch data id
    */
   int getCPatchDataId() const;

   /*!
    * @brief Get C as a constant value.
    *
    * Error if C is not represented by a constant.
    *
    * @return C's constant value
    */
   double getCConstant() const;

   //@}


private:

   /*!
    * @brief Object name.
    */
   std::string d_object_name;

   int d_D_id;
   double d_D_constant;

   bool d_C_zero;
   int d_C_id;
   double d_C_constant;

};

} // namespace SAMRAI
}

#ifndef DEBUG_NO_INLINE
#include "PoissonSpecifications.I"
#endif

#endif
