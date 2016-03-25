//
// File:        solv::PVodeTrioAbstractVector.h
// Package:     SAMRAI solvers package
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to C++ vector kernel operations for PVodeTrio package.
//

#ifndef included_solv_PVodeTrioAbstractVector
#define included_solv_PVodeTrioAbstractVector

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#if HAVE_KINSOL || HAVE_PVODE


namespace SAMRAI {
   namespace solv {


/**
 * Class solv::PVodeTrioAbstractVector is an abstract base class that declares 
 * operations provided by any {\tt C++} class that may be used as the
 * vector kernel by the PVodeTrio nonlinear solver package.  PVodeTrio allows 
 * arbitrarily defined vectors to be used within it as long as the proper
 * collection of operations are provided.  The intent of this base class 
 * to provide the interface for one's own vector kernel.  One implements
 * a subclass that provides the functions declared herein as pure virtual 
 * and which provides necessary vector data structures.  Note that the 
 * virtual members of this class are all protected.  They should not be used 
 * outside of a subclass of this class.
 *
 * PVodeTrio was developed in the Center for Applied Scientific Computing (CASC)
 * at Lawrence Livermore National Laboratory (LLNL).  For more information
 * about PVodeTrio, see A.G. Taylor and A.C. Hindmarsh, "User documentation for
 * PVodeTrio, a nonlinear solver for sequential and parallel computers",
 * UCRL-ID-131185, Lawrence Livermore National Laboratory, 1998.
 *
 * Important notes:
 * 


 * - \b (1) The user-supplied vector subclass should only inherit from
 *            this base class. It MUST NOT employ multiple inheritance so that
 *            problems with casting from a base class pointer to a subclass
 *            pointer and the use of virtual functions works properly.
 * - \b (2) The user-supplied subclass that implements the vector data
 *            structures and operations is responsible for all parallelism,
 *            I/O, etc. associated with the vector objects.  PVodeTrio is 
 *            implemented using a SPMD programming model.  It has no knowlege
 *            of the structure of the vector data, nor the implementations
 *            of the individual vector routines.
 * - \b (3) We assume the vector data is {\tt double}, which is the
 *            default for PVodeTrio.
 * 


 *
 * @see solv::PVodeTrioSolver
 */

class PVodeTrioAbstractVector
{
public:
   /**
    * Uninteresting constructor and destructor for solv::PVodeTrioAbstractVector.
    */
   PVodeTrioAbstractVector();
   virtual ~PVodeTrioAbstractVector(); 

   /**
    * Clone the vector structure and allocate storage for the vector
    * data.  Then, return a pointer to the new vector instance.  Note that 
    * the new vector object must be distinct from the original.  This
    * function is distinct from the vector constructor since it will
    * be called from within PVodeTrio to allocate vectors during the nonlinear
    * solution process.  The original solution vector must be setup by the 
    * user's application code.
    */
   virtual PVodeTrioAbstractVector* makeNewVector() = 0;

   /**
    * Destroy vector structure and its storage. This function is distinct 
    * from the destructor since it will be called from within PVodeTrio to 
    * deallocate vectors during the nonlinear solution process.
    */
   virtual void freeVector() = 0;

   /**
    * Initialize all entries of this vector object to scalar \f$c\f$.
    */
   virtual void setToScalar(const double c) = 0;

   /**
    * Set this vector object to scalar \f$c x\f$, where \f$c\f$ is a scalar and
    * x is another vector.
    */
   virtual void scaleVector(const PVodeTrioAbstractVector* x,
                            const double c) = 0;

   /**
    * Set this vector object to \f$a x + b y\f$, where \f$a, b\f$ are scalars and
    * \f$x, y\f$ are vectors. 
    */ 
   virtual void setLinearSum(const double a,
                             const PVodeTrioAbstractVector* x, 
                             const double b,
                             const PVodeTrioAbstractVector* y) = 0;

   /**
    * Set each entry of this vector: \f$v_i = x_i y_i\f$, where \f$x_i, y_i\f$ are
    * entries in vectors \f$x\f$ and \f$y\f$.
    */
   virtual void pointwiseMultiply(const PVodeTrioAbstractVector* x,
                                  const PVodeTrioAbstractVector* y) = 0;

   /**
    * Set each entry of this vector: \f$v_i = x_i / y_i\f$, where \f$x_i, y_i\f$ are
    * entries in vectors \f$x\f$ and \f$y\f$.  Based on the PVodeTrio vector 
    * implementation, it is not necessary to check for division by zero.
    */
   virtual void pointwiseDivide(const PVodeTrioAbstractVector* x,
                                const PVodeTrioAbstractVector* y) = 0;

   /**
    * Set each entry of this vector to the absolute value of the 
    * corresponding entry in vector \f$x\f$.
    */
   virtual void setAbs(const PVodeTrioAbstractVector* x) = 0;

   /**
    * Set each entry of this vector: \f$v_i =  1 / x_i\f$, where \f$x_i\f$ is an entry
    * entry in vector \f$x\f$.  Based on the PVodeTrio vector implementation,
    * it is not necessary to no check for division by zero.
    */
   virtual void pointwiseReciprocal(const PVodeTrioAbstractVector* x) = 0;

   /**
    * Set each entry of this vector to the corresponding entry in vector \f$x\f$
    * plus the scalar \f$b\f$.
    */
   virtual void addScalar(const PVodeTrioAbstractVector* x,
                          const double b) = 0;

   /**
    * Return the dot product of this vector and the argument vector \f$x\f$.
    */
   virtual double dotWith(const PVodeTrioAbstractVector* x) const = 0;

   /**
    * Return the max norm of this vector.
    */
   virtual double maxNorm() const = 0;

   /**
    * Return the \f$L_1\f$ norm of this vector.
    */
   virtual double L1Norm() const = 0;

   /**
    * Return the weighted-\f$L_2\f$ norm of this vector using the vector \f$x\f$ 
    * as the weighting vector.
    */
   virtual double weightedL2Norm(const PVodeTrioAbstractVector* x) const = 0;

   /**
    * Return the weighted root mean squared norm of this vector using 
    * the vector \f$x\f$ as the weighting vector.
    */
   virtual double weightedRMSNorm(
      const PVodeTrioAbstractVector* x) const = 0;

   /**
    * Return the minimum entry of this vector.
    */
   virtual double vecMin() const = 0;

   /**
    * Return \f$0\f$ if \f$x_i \neq 0\f$ and \f$x_i z_i \leq 0\f$, for some \f$i\f$.
    * Here \f$z_i\f$ is an element of this vector. Otherwise, return \f$1\f$.
    */
   virtual int constrProdPos(const PVodeTrioAbstractVector* x) const = 0;

   /**
    * Set each entry in this vector based on the vector \f$x\f$ as follows:
    * if \f$\mid x_i \mid \geq c\f$, then \f$v_i = 1\f$, else \f$v_i = 0\f$.
    */
   virtual void compareToScalar(const PVodeTrioAbstractVector* x,
                                const double c) = 0;

   /**
    * Set each entry of this vector: \f$v_i =  1 / x_i\f$, where \f$x_i\f$ is an
    * entry entry in the vector \f$x\f$, unless \f$x_i = 0\f$.  If \f$x_i = 0\f$,
    * then return \f$0\f$.  Otherwise, \f$1\f$ is returned.
    */
   virtual int testReciprocal(const PVodeTrioAbstractVector* x) = 0;

   /**
    * Print the vector data to the output stream used by the subclass 
    * print routine. 
    */
   virtual void printVector() const = 0;
};


}
}
#endif

#endif
