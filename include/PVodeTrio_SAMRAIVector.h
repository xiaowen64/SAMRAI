//
// File:        PVodeTrio_SAMRAIVector.h
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: "Glue code" between PVodeTrio vector interface and SAMRAI vectors.
//

#ifndef included_solv_PVodeTrio_SAMRAIVector
#define included_solv_PVodeTrio_SAMRAIVector

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT     
*  KINSOL -or- PVODE
************************************************************************
*/
#if HAVE_KINSOL || HAVE_PVODE

#ifndef included_solv_PVodeTrioAbstractVector
#include "PVodeTrioAbstractVector.h"
#endif
#ifndef included_solv_SAMRAIVectorReal
#include "SAMRAIVectorReal.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace solv {

/**
 * Class PVodeTrio_SAMRAIVector<DIM> wraps a real-valued SAMRAI vector 
 * (see SAMRAIVectorReal class) object so that it may be used with 
 * the PVodeTrio solver packages.  This class is derived from the 
 * abstract base class PVodeTrioAbstractVector, which defines a {\tt C++} 
 * interface for PVodeTrio vectors.  It also maintains a pointer to a SAMRAI 
 * vector object.  A SAMRAI vector is defined as a collection of patch data 
 * components living on some subset of levels in a structured AMR mesh 
 * hierarchy.
 *
 * Observe that there are only three public member functions in this class
 * They are used to create and destroy PVodeTrio vector objects (i.e., 
 * "N_Vector"s), and to obtain the SAMRAI vector associated with the PVodeTrio 
 * vector.  In particular, note that the constructor and destructor of this 
 * class are protected members.  The construction and destruction of instances 
 * of this class may occur only through the static member functions that
 * create and destroy PVodeTrio vector objects.
 *
 * Finally, we remark that this class provides vectors of type {\tt double},
 * which is the default for PVodeTrio.
 * 
 * @see solv::PVodeTrioAbstractVector
 * @see solv::SAMRAIVectorReal 
 */

template<int DIM> class PVodeTrio_SAMRAIVector : public PVodeTrioAbstractVector
{
public:
   /**
    * Create and return a new PVodeTrio vector object.  The SAMRAI vector
    * object is wrapped so that it may be manipulated within PVodeTrio
    * as an N_Vector (which is typedef'd to PVodeTrioAbstractVector* 
    * in the abstract PVodeTrio vector interface).  It is important to note 
    * that this function does not allocate storage for the vector data.  
    * Data must be allocated through the SAMRAI vector object directly.
    * For output of the data through "N_VPrint" calls, the output stream 
    * to which the SAMRAI vector object writes will be used.
    */
   static PVodeTrioAbstractVector* createPVodeTrioVector(
      tbox::Pointer< SAMRAIVectorReal<DIM,double> > samrai_vec);

   /**
    * Destroy a given PVodeTrio vector object. It is important to note that
    * this function does not deallocate storage for the vector data.
    * Vector data must be deallocated hrough the SAMRAI vector object.
    */
   static void destroyPVodeTrioVector(
      PVodeTrioAbstractVector* pvode_trio_vec);

   /**
    * Return pointer to the SAMRAI vector object associated with the
    * given PVodeTrio vector.
    */
   static tbox::Pointer< SAMRAIVectorReal<DIM,double> > 
   getSAMRAIVector(PVodeTrioAbstractVector* pvode_trio_vec);

protected:
   /*
    * Constructor for PVodeTrio_SAMRAIVector<DIM>.
    */
   PVodeTrio_SAMRAIVector(
      tbox::Pointer< SAMRAIVectorReal<DIM,double> > samrai_vector);

   /*
    * Virtual destructor for PVodeTrio_SAMRAIVector<DIM>.
    */
   virtual ~PVodeTrio_SAMRAIVector<DIM>();

private:
   /*
    * Return SAMRAI vector owned by this PVodeTrio_SAMRAIVector<DIM> object.
    */
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > getSAMRAIVector();

   /*
    * The makeNewVector() function clones the vector structure and allocate 
    * storage for the new vector data.  Then, a pointer to the new vector 
    * instance is returned. This function is distinct from the constructor 
    * since it will be called from within PVodeTrio to allocate vectors used
    * during the nonlinear solution process.
    */
   PVodeTrioAbstractVector* makeNewVector();

   /*
    * Destroy vector structure and its storage. This function is distinct 
    * from the destructor since it will be called from within PVodeTrio to 
    * deallocate vectors during the nonlinear solution process.
    */
   virtual void freeVector();

   /*
    * Initialize all entries of this vector object to scalar \f$c\f$.
    */
   void setToScalar(const double c);

   /*
    * Set this vector object to scalar \f$c x\f$, where \f$c\f$ is a scalar and
    * \f$x\f$ is another vector.
    */
   void scaleVector(const PVodeTrioAbstractVector* x,
                    const double c);

   /*
    * Set this vector object to \f$a x + b y\f$, where \f$a, b\f$ are scalars and
    * \f$x, y\f$ are vectors. 
    */ 
   void setLinearSum(const double a, 
                     const PVodeTrioAbstractVector* x,
                     const double b, 
                     const PVodeTrioAbstractVector* y);

   /*
    * Set each entry of this vector: \f$v_i = x_i y_i\f$, where \f$x_i, y_i\f$ are
    * entries in vectors \f$x\f$ and \f$y\f$.
    */
   void pointwiseMultiply(const PVodeTrioAbstractVector* x,
                          const PVodeTrioAbstractVector* y);

   /*
    * Set each entry of this vector: \f$v_i = \frac{x_i}{y_i}\f$, where 
    * \f$x_i, y_i\f$ are entries in vectors \f$x\f$ and \f$y\f$.  Based on the PVodeTrio 
    * vector implementation, it is not necessary to no check for division by 
    * zero.
    */
   void pointwiseDivide(const PVodeTrioAbstractVector* x,
                        const PVodeTrioAbstractVector* y);

   /*
    * Set each entry of this vector to the absolute value of the 
    * corresponding entry in vector \f$x\f$.
    */
   void setAbs(const PVodeTrioAbstractVector* x);

   /*
    * Set each entry of this vector: \f$v_i = \frac{1}{x_i}\f$, where \f$x_i\f$ is 
    * an entry in vector \f$x\f$.  Based on the PVodeTrio vector implementation,
    * it is not necessary to no check for division by zero.
    */
   void pointwiseReciprocal(const PVodeTrioAbstractVector* x);

   /*
    * Set each entry of this vector: \f$v_i = x_i + b\f$, where \f$x_i\f$ is an entry 
    * in the vector \f$x\f$ and \f$b\f$ is a scalar.
    */
   void addScalar(const PVodeTrioAbstractVector* x,
                  const double b);

   /*
    * Return the dot product of this vector and the vector \f$x\f$.
    */
   double dotWith(const PVodeTrioAbstractVector* x) const;

   /*
    * Return the max norm of this vector: 
    * \f${\| v \|}_{\max} = \max_{i} (\mid v_i \mid)\f$.
    */
   double maxNorm() const;

   /*
    * Return the \f$L_1\f$ norm of this vector: 
    * \f${\| v \|}_{L_1} = \sum_{i} (\mid v_i \mid)\f$ if no control volumes
    * are defined.  Otherwise, 
    * \f${\| v \|}_{L_1} = \sum_{i} (\mid v_i \mid * cvol_i )\f$.
    */
   double L1Norm() const;

   /*
    * Return the weighted \f$L_2\f$ norm of this vector using the vector 
    * \f$x\f$ as the weighting vector: 
    * \f${\| v \|}_{WL2(x)} = \sqrt{ \sum_i( (x_i * v_i)^2 ) )}\f$ if no
    * control volumes are defined.  Otherwise, 
    * \f${\| v \|}_{WL2(x)} = \sqrt{ \sum_i( (x_i * v_i)^2 cvol_i ) }\f$.
    */
   double weightedL2Norm(const PVodeTrioAbstractVector* x) const;

   /*
    * Return the weighted root mean squared norm of this vector using 
    * the vector \f$x\f$ as the weighting vector.   If control volumes are 
    * not defined for the vector entries, the norm corresponds to the 
    * weighted \f$L_2\f$-norm divided by the square root of the number of 
    * vector entries.  Otherwise, the norm corresponds to the weighted 
    * \f$L_2\f$-norm divided by the square root of the sum of the control volumes.
    */
   double weightedRMSNorm(const PVodeTrioAbstractVector* x) const;

   /*
    * Return the minimum entry of this vector.
    */
   double vecMin() const;

   /*
    * Return \f$0\f$ if \f$x_i \neq 0\f$ and \f$x_i v_i \leq 0\f$, for some \f$i\f$.
    * Otherwise, return \f$1\f$.
    */
   int constrProdPos(const PVodeTrioAbstractVector* x) const;

   /*
    * Set each entry in this vector based on the vector \f$x\f$ as follows:
    * if \f$\mid x_i \mid \geq c\f$, then \f$v_i = 1\f$, else \f$v_i = 0\f$.
    */
   void compareToScalar(const PVodeTrioAbstractVector* x,
                        const double c);

   /*
    * Set each entry of this vector: \f$v_i =  \frac{1}{x_i}\f$, where \f$x_i\f$ 
    * is an entry entry in the vector \f$x\f$, unless \f$x_i = 0\f$.  If \f$x_i = 0\f$, 
    * then return \f$0\f$.  Otherwise, \f$1\f$ is returned.
    */
   int testReciprocal(const PVodeTrioAbstractVector* x);

   /*
    * Print the vector to the output stream used by the SAMRAI vector class.
    */
   void printVector() const;

   /*
    * Vector data is maintained in SAMRAI vector structure.
    */
   tbox::Pointer< SAMRAIVectorReal<DIM,double> > d_samrai_vector;
};

}
}
#ifndef DEBUG_NO_INLINE
#include "PVodeTrio_SAMRAIVector.I"
#endif
#endif

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PVodeTrio_SAMRAIVector.C"
#endif
