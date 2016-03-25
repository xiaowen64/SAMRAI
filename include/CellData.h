//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/cell/CellData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated cell centered patch data type
//

#ifndef included_pdat_CellData
#define included_pdat_CellData

#include "SAMRAI_config.h"
#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#include "PatchData.h"
#include "ArrayData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "tbox/Arena.h"
#include "tbox/Complex.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class CellData<DIM> provides an implementation for data defined 
 * at cell centers on AMR patches.  It is derived from the 
 * hier::PatchData interface common to all SAMRAI patch data types.  
 * Given a CELL-centered AMR index space box, a cell data object represents 
 * data of some template TYPE and depth at the centers of the cells in the box.
 * Here, depth indicates the number of data values at each cell index location.
 * The CellGeometry class provides the translation between the standard SAMRAI 
 * cell-centered AMR index space and cell-centered data.
 *
 * A cell-centerd data array is stored in (i,...,k,d) order, where 
 * i,...,k are spatial indices and d indicates the depth at that location.
 * Memory allocation is in column-major ordering (e.g., Fortran style)
 * so that the leftmost index runs fastest in memory.  For example, a 
 * three-dimensional cell data object defined over a box 
 * [l0:u0,l1:u1,l2:u2] holds a data array dimensioned as 
 * \verbatim

     [ l0 : u0 ,
       l1 : u1 ,
       l2 : u2 , d ]

 * \endverbatim
 * Other spatial dimensions are represented similarly.
 *
 * The data type TYPE must define a default constructor (i.e., taking no
 * arguments) and also the copy assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::CellDataFactory
 * @see pdat::CellIndex
 * @see pdat::CellIterator
 * @see pdat::CellGeometry
 */

template<int DIM, class TYPE>
class CellData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent cell-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            cell data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the node data will be allocated.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box, 
                               int depth, 
                               const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief The constructor for an cell data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            cell data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *              of the ghost cell region around the box over which
    *              the node data will be allocated.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   CellData(const hier::Box<DIM>& box,
            int depth,
            const hier::IntVector<DIM>& ghosts,
            tbox::Pointer<tbox::Arena> pool = tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief The virtual destructor for a cell data object.
    */
   virtual ~CellData<DIM,TYPE>();

   /*!
    * @brief Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int getDepth() const;
 
   /*!
    * @brief Get a pointer to the beginning of a depth 
    * component of the cell centered array.
    */
   TYPE* getPointer(int depth = 0);
 
   /*!
    * @brief Get a const pointer to the beginning of a depth 
    * component of the cell centered array.
    */
   const TYPE* getPointer(int depth = 0) const;

   /*!
    * @brief Return reference to cell data entry corresponding
    * to a given cell index and depth.
    */
   TYPE& operator()(const CellIndex<DIM>& i, int depth = 0);

   /*!
    * @brief Return a const reference to cell data entry corresponding
    * to a given cell index and depth.
    */
   const TYPE& operator()(const CellIndex<DIM>& i, int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for
    * the cell centered data object.
    */
   ArrayData<DIM,TYPE>& getArrayData();
 
   /*!
    * @brief Return a const reference to the array data object for
    * the cell centered data object.
    */
   const ArrayData<DIM,TYPE>& getArrayData() const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * CellData of the same DIM and TYPE.  If not, then an unrecoverable 
    * error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * CellData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be CellData of the same DIM and TYPE
    * and the overlap must be a CellOverlap of the same DIM.  
    * If not, then an unrecoverable error results.
    */
   virtual void copy(const hier::PatchData<DIM>& src, 
                     const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be CellData of the same DIM and TYPE
    * and the overlap must be a CellOverlap of the same DIM.  
    * If not, then an unrecoverable error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst, 
                      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given CELL-centered AMR index box.
    */
   void copyOnBox(const CellData<DIM,TYPE>& src, 
                  const hier::Box<DIM>& box);

   /*!
    * @brief Fast copy (i.e., source and this cell data objects are
    * defined over the same box) to this destination cell data object  
    * from the given source cell data object at the specified depths.
    */
   virtual void copyDepth(int dst_depth,
			  const CellData<DIM,TYPE>& src,
			  int src_depth);

   /*!
    * @brief Return true if the patch data object can estimate the
    * stream size required to fit its data using only index
    * space information (i.e., a box).
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual bool canEstimateStreamSizeFromBox() const;

   /*!
    * @brief Return the number of bytes needed to stream the data
    * in this patch data object lying in the specified box overlap
    * region.
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual int getDataStreamSize(const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region.  The overlap must be a 
    * CellOverlap of the same DIM.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object
    * over the specified box overlap region.  The overlap must be a  
    * CellOverlap of the same DIM.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Fill all values at depth d with the value t.
    */
   void fill(const TYPE& t, 
             int d = 0);

   /*!
    * @brief Fill all values at depth d within the box with the value t.
    */
   void fill(const TYPE& t, 
             const hier::Box<DIM>& box, 
             int d = 0);

   /*!
    * @brief Fill all depth components with value t.
    */
   void fillAll(const TYPE& t);

   /*!
    * @brief Fill all depth components within the box with value t.
    */
   void fillAll(const TYPE& t, 
                const hier::Box<DIM>& box);

   /*!
    * @brief Print all cell data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to cell index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void print(const hier::Box<DIM>& box, 
              std::ostream& os = tbox::plog, 
              int prec = 12) const;

   /*!
    * @brief Print all cell data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to cell index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void print(const hier::Box<DIM>& box, 
              int depth, 
              std::ostream& os = tbox::plog, 
              int prec = 12) const;

   /*!
    * Check that class version and restart file version are equal.  
    * If so, read data members from the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase( 
        tbox::Pointer<tbox::Database> database);

   /*!
    * Write out the class version number and other data members to 
    * the database. 
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void putSpecializedToDatabase( 
        tbox::Pointer<tbox::Database> database);

   /*!
    * The cell iterator iterates over the elements of a cell
    * centered box geometry.  This typedef is a convenience 
    * for using the CellIterator<DIM> class.
    */
   typedef CellIterator<DIM> Iterator;

private:
   CellData(const CellData<DIM,TYPE>&);	// not implemented
   void operator=(const CellData<DIM,TYPE>&);		// not implemented

   int d_depth;
   ArrayData<DIM,TYPE> d_data;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "CellData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellData.C"
#endif
