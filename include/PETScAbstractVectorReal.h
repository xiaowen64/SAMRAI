//
// File:        PETScAbstractVectorReal.h
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface to C++ vector implementation for PETSc package.
//

#ifndef included_solv_PETScAbstractVectorReal
#define included_solv_PETScAbstractVectorReal

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT PETSC
************************************************************************
*/
#ifdef HAVE_PETSC

#ifdef REQUIRES_CMATH
// SGS include for blue xlc
#include <cmath>
#endif

#ifndef included_petsc_vec
#define included_petsc_vec
extern "C" {
#include "petscvec.h"
}
#endif

namespace SAMRAI {
    namespace solv {

/**
 * Class PETScAbstractVectorReal serves as an abstract base class for a
 * {\tt C++} vector class that can be used with the PETSc solver framework.
 * Specifically, this class provides an interface for real-valued PETSc
 * vectors (i.e., where the data is either float or double).  PETSc allows 
 * the use of user-defined vectors.  Thus, the intent of this base class is 
 * that one may provide his/her own vector implementation in a subclass of 
 * this base class that provides the the necessary vector data structures
 * and implements the pure virtual functions.  This class declares private 
 * static member functions for linkage with the vector object function calls 
 * understood by PETSc.  Each of the static member functions calls an 
 * associated function in a subclass via the virtual function mechanism.  
 * Note that the virtual members of this class are all protected.  They should 
 * not be used outside of a subclass of this class.  The data member of the 
 * PETSc vector object is set to an instantiation of the user-supplied vector 
 * class, when an object of this class is constructed. 
 *
 * PETSc was developed in the Mathematics and Computer Science Division at 
 * Argonne National Laboratory (ANL).  For more information about PETSc, 
 * see {\tt http://www-fp.mcs.anl.gov/petsc/}.
 *
 * Important notes: 
 * 


 * - @b (1) The user-supplied vector subclass should only inherit from 
 *            this base class. It MUST NOT employ multiple inheritance so that
 *            problems with casting from a base class pointer to a subclass
 *            pointer and the use of virtual functions works properly. 
 * - @b (2) The user-supplied subclass that implements the vector data 
 *            structures and operations is responsible for all parallelism, 
 *            I/O, etc. associated with the vector objects.  PETSc only sees 
 *            pointers to what it believes to be sequential vector objects 
 *            associated with the local process only.  It has no knowledge
 *            of the structure of the vector data, nor the implementations
 *            of the individual vector routines.
 * - @b (3) Several of the operations defined in the PETSc {\tt _VecOps}
 *            structure are left unimplemented in this class.  They will
 *            print an error message and throw an unrecoverable exeception
 *            if called which causes the program to abort. 
 * - @b (4) By default, PETSc typdefs "Scalar" to {\tt double}.  PETSc
 *            must be recompiled to use {\tt float} data.  Also, PETSc
 *            support complex vector data.  A complex vector interface class 
 *            similar to this class may be provided in the future if the 
 *            need arises.
 * 


 */

template <class TYPE>
class PETScAbstractVectorReal
{
protected:
   /**
    * Constructor PETScAbstractVectorReal class that provides a wrapper
    * for a SAMRAI vector so that it can be manipulated within PETSc.  The
    * constructor allocates the PETSc vector and sets its data structures and
    * member functions so that it can operate on the SAMRAI vector.
    */
   PETScAbstractVectorReal(bool vector_created_via_duplicate);

   /**
    * Destructor for PETScAbstractVectorReal class destroys the PETSc
    * vector created by the constructor.
    */ 
   virtual ~PETScAbstractVectorReal(); 

   /**
    * Return PETSc "Vec" object for this PETScAbstractVectorReal object.
    */
   Vec getPETScVector();

   /**
    * Clone the vector structure and allocate storage for the vector
    * data.  Then, return a pointer to the new vector instance.  This
    * function is distinct from the vector constructor since it is called
    * from within PETSc (via the duplicateVec(), duplicateVecs() functions)
    * to allocate new vector objects.
    */
   virtual PETScAbstractVectorReal<TYPE>* makeNewVector() = 0;

   /**
    * Destroy vector structure and its storage. This function is distinct
    * from the destructor since it will be called from within PETSc to
    * deallocate vectors (via the destroyVec(), destroyVecs() functions). 
    */
   virtual void freeVector() = 0;

   /**
    * View all vector data.  Note that the user-supplied vector must 
    * choose how to view the vector and its data; e.g., print to file, 
    * dump to standard out, etc.
    */
   virtual void viewVector() const = 0;

   /**
    * Return @f$ (x,y) = \sum_i ( x_i * conj(y_i) ) @f$ , where @f$ x @f$  is this vector.
    * Note that for real vectors, this is the same as TdotWith().
    * If local_only is true, the operation is limited to parts owned by the
    * local process.
    */
   virtual TYPE dotWith(const PETScAbstractVectorReal<TYPE>* y,
			bool local_only=false) const = 0;

   /**
    * Limited to local data only,
    * return @f$ (x,y) = \sum_i ( x_i * (y_i) ) @f$ , where @f$ x @f$  is this vector.
    * Note that for real vectors, this is the same as dotWith().
    * If local_only is true, the operation is limited to parts owned by the
    * local process.
    */
   virtual TYPE TdotWith(const PETScAbstractVectorReal<TYPE>* y,
			bool local_only=false) const = 0;

   /**
    * Return @f$ L_1 @f$ -norm of this vector.
    *
    * @param local_only Flag to get result for local data only.
    */
   virtual double L1Norm(bool local_only=false) const = 0;

   /**
    * Return @f$ L_2 @f$ -norm of this vector.
    *
    * @param local_only Flag to get result for local data only.
    */
   virtual double L2Norm(bool local_only=false) const = 0;

   /**
    * Return @f$ L_{\infty} @f$ -norm of this vector.
    *
    * @param local_only Flag to get result for local data only.
    */
   virtual double maxNorm(bool local_only=false) const = 0;

   /**
    * Multiply each entry of this vector by given scalar.
    */
   virtual void scaleVector(const TYPE alpha) = 0;

   /**
    * Copy source vector data to this vector.
    */
   virtual void copyVector(const PETScAbstractVectorReal<TYPE>* v_src) = 0;

   /**
    * Set each entry of this vector to given scalar.
    */
   virtual void setToScalar(const TYPE alpha) = 0;

   /**
    * Swap data between this vector and argument vector.
    */
   virtual void swapWith(PETScAbstractVectorReal<TYPE>* v_other) = 0; 

   /**
    * Set @f$ y = \alpha x + y @f$ , where @f$ y @f$  is this vector.
    */
   virtual void setAXPY(const TYPE alpha,
                        const PETScAbstractVectorReal<TYPE>* x) = 0;

   /**
    * Set @f$ y = \alpha x + @beta y @f$ , where @f$ y @f$  is this vector.
    */
   virtual void setAXPBY(const TYPE alpha,
                         const PETScAbstractVectorReal<TYPE>* x,
                         const TYPE beta) = 0;

   /**
    * Set @f$ w = \alpha x + y @f$ , where @f$ w @f$  is this vector.
    */
   virtual void setWAXPY(const TYPE alpha,
                         const PETScAbstractVectorReal<TYPE>* x,
                         const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Set @f$ w_i = x_i y_i @f$ , where @f$ w_i @f$  is @f$ i @f$ -th entry of this vector.
    */
   virtual
   void pointwiseMultiply(const PETScAbstractVectorReal<TYPE>* x,
                          const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Set @f$ w_i = x_i / y_i @f$ , where @f$ w_i @f$  is @f$ i @f$ -th entry of this vector.
    */
   virtual
   void pointwiseDivide(const PETScAbstractVectorReal<TYPE>* x,
                        const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Compute @f$ max_i = abs(w_i / y_i) @f$ ,
    where @f$ w_i @f$  is @f$ i @f$ -th entry of this vector.
    */
   virtual
   double maxPointwiseDivide(const PETScAbstractVectorReal<TYPE>* y) = 0;

   /**
    * Find maximum vector entry and vector index at which maximum occurs.
    */
   virtual void vecMax(int& i, double& max) const = 0;

   /**
    * Find minimum vector entry and vector index at which minimum occurs.
    */
   virtual void vecMin(int& i, double& min) const = 0;

   /**
    * Set vector entries to random values.  Note that PETSc uses the
    * drand48() function and computes random value as width*drand48()+low.
    */
   virtual void setRandomValues(const TYPE width, const TYPE low) = 0;

   /**
    * Set argument to vector data in contiguous array (local to processor).
    */
   virtual void getDataArray(TYPE* array) = 0;

   /**
    * Return total length of vector. 
    */
   virtual int getDataSize() const = 0;

   /**
    * Return length of vector (local to processor).
    */
   virtual int getLocalDataSize() const = 0;

private:
   /*
    * PETSc vector object corresponding to this 
    * PETScAbstractVectorReal object.
    */
   Vec d_petsc_vector; 
   
   bool d_vector_created_via_duplicate;

   /*
    * Static member functions for linkage with PETSc solver package routines.
    * Essentially, these functions match those in the PETSc _VecOps structure.
    * Note that these operations are actually implemented in a subclass of 
    * this base class using the virtual function mechanism.
    */

   /* Free vector structure and associated data. */
   static int destroyVec(Vec v);

   /* Print vector data. */
   static int viewVec(Vec x, PetscViewer view);

   /* Duplicate vector structure and allocate data storage for new vector. */ 
   static int duplicateVec(Vec v_in, Vec* v_new);

   /* Duplicate array of vectors and allocate data storage for new vectors. */ 
   static int duplicateVecs(Vec v_in, int n, Vec** varr_new);

   /* Free each vector structure and its data in vector array. */
   static int destroyVecs(const Vec* v_arr, int n);

   /* Compute dot product: *dp = (x,y) = sum( x_i * conj(y_i) ) */
   static int dotProduct(Vec x, Vec y, TYPE* dp);

   /* Compute dot product: *dp = (x,y)_T = sum( x_i * y_i ) */
   static int dotProductT(Vec x, Vec y, TYPE* dp);

   /* Compute dot product: dp[j] = (x,y_j) = sum( x_i * conj((y_j)_i) ) */
   static int dotProductM(int n, Vec x, const Vec* y, TYPE* dp);

   /* Compute dot product: dp[j] = (x,y_j)_T = sum( x_i * (y_j)_i ) */
   static int dotProductMT(int n, Vec x, const Vec* y, TYPE* dp);

   /*
      Compute requested norm, where PETSc defines the following enumerated
      type (in vec.h):
            typedef enum {NORM_1=1,
                          NORM_2=2,
                          NORM_FROBENIUS=3,
                          NORM_INFINITY=4,
                          NORM_1_AND_2=5} NormType;
      If norm type is not NORM_1, NORM_2, NORM_INFINITY, or NORM_1_AND_2
      an unrecoverable exception will be thrown and program will abort. 
   */ 
   static int vecNorm(Vec x, NormType n_type, double* norm);

   /* Scale vector entries: x = alpha * x   */
   static int scaleVec(const TYPE* alpha, Vec x);

   /* Copy source vector to destination: v_dst = v_src */
   static int copyVec(Vec v_src, Vec v_dst);

   /* Set vector entries to scalar: x = alpha  */
   static int setVec(const TYPE* alpha, Vec x);

   /* Exchange vectors x and y. */
   static int swapVecs(Vec x, Vec y);

   /* Set y = alpha * x + y */
   static int computeAXPY(const TYPE* alpha, Vec x, Vec y);

   /* Set y = alpha * x + beta * y */
   static int computeAXPBY(const TYPE* alpha, const TYPE* beta, Vec x, Vec y); 

   /* Set x = x + SUM_j (alpha[j] + y[j]) */
   static int computeMAXPY(int n, const TYPE* alpha, Vec x, Vec* y); 

   /* Set y = x + alpha * y */
   static int computeAYPX(const TYPE* alpha, Vec x, Vec y);

   /* Set w = alpha * x + y */
   static int computeWAXPY(const TYPE* alpha, Vec x, Vec y, Vec w);

   /* Set w_i = x_i * y_i */
   static int pointwiseMultVecs(Vec x, Vec y, Vec w);

   /* Set w_i = x_i / y_i */
   static int pointwiseDivideVecs(Vec x, Vec y, Vec w);

   /* Set max = max of abs(x_i / y_i) */
   static int maxPointwiseDivideVecs(Vec x, Vec y, PetscReal *max);

   /* Compute *max = max_i(x_i) and  *i = index of max elt. */
   static int vectorMax(Vec x, int* i, double* max);

   /* Compute *min = min_i(x_i) and *i = index of min elt. */
   static int vectorMin(Vec x, int* i, double* min);

   /* Set vector entries to random values.  Note: PETSc uses drand48().  */
   static int setVecRandom(PetscRandom r, Vec x);

   /* Return *array = vector data in contiguous array (local to processor). */
   static int getVecArray(Vec x, TYPE** array);

   /* Return *size = total number of vector entries. */
   static int getVecSize(Vec x, int* size);

   /* Return *size = number of local vector entries. */
   static int getLocalVecSize(Vec x, int* size);

   /*
    *********************************************************************
    * The remaining functions are not implemented and will result in an *  
    * unrecoverable exception being thrown and program abort if called. *
    *********************************************************************
    */
 
   /*
      For each vector elt index in indices (size = ni),
      set x(index) = vals(index) if m == INSERT_VALUES.
      If m != INSERT_VALUES, set x(index) += vals(index).
      Note: PETSc routine checks for index out of range.
            Also, enum InsertMode mode type is defined in vec.h
            (has several modes defined).
   */
   static int setVecValues(Vec x, int ni, const int* indices, const TYPE* vals,
                           InsertMode mode);

   static int beginVecAssembly(Vec x);

   static int endVecAssembly(Vec x);

   /* Restore vector data array. */
   static int restoreVecArray(Vec x, TYPE** array);

   static int setVecOption(Vec x, VecOption option);

   static int setVecValuesBlocked(Vec x, int n, const int* nb, 
                                  const TYPE* vals,
                                  InsertMode mode); 

   static int placeArray(Vec x, const TYPE *vals);

   static int replaceArray(Vec x, const TYPE *vals);

   static int dotProductLocal(Vec x, Vec y, TYPE* dp);

   static int dotProductLocalT(Vec x, Vec y, TYPE* dp);

   static int vecNormLocal(Vec x, NormType n_type, double* norm);

   static int loadIntoVector(PetscViewer view, Vec x);

   static int reciprocal(Vec x);

   static int viewNative(Vec x, PetscViewer view);

   static int conjugate(Vec x);

   static int setLocalToGlobalMapping(Vec x, ISLocalToGlobalMapping mapping);

   static int setValuesLocal(Vec x, int count, const int *indices,
                             const TYPE *vals, InsertMode mode);

   static int resetArray(Vec x);

   static int setFromOptions(Vec x);
};

}
}
#endif
#endif
