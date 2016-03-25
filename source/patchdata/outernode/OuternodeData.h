//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/outernode/OuternodeData.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Release:	$Name$
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Templated outernode centered patch data type
//

#ifndef included_pdat_OuternodeData
#define included_pdat_OuternodeData

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_pdat_ArrayData
#include "ArrayData.h"
#endif
#ifndef included_pdat_NodeData
#include "NodeData.h"
#endif
#ifndef included_pdat_NodeIndex
#include "NodeIndex.h"
#endif
#ifndef included_pdat_NodeOverlap
#include "NodeOverlap.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class OuternodeData<DIM> provides an implementation for data defined
 * at cell nodes on the boundaries of AMR patches.  It is derived from the
 * hier::PatchData interface common to all SAMRAI patch data types.  Given
 * a CELL-centered AMR index space box, an outernode data object represents
 * data of some template TYPE and depth on the cell nodes on the boundary
 * of the box.  Here, depth indicates the number of data values at each node
 * index location.  The OuternodnodeGeometry class provides the translation
 * between the standard SAMRAI cell-centered AMR index space and
 * outernode-centered data.
 *
 * Outernode data is stored in 2*DIM arrays, each of which contains data
 * associated with node indices on an upper or lower box face in some
 * corrdinate direction.  The data layout in the outernode data arrays matches
 * the corresponding array sections provided by the node data implementation.
 * Where a node index falls on more than one box face (patch boundary edges and 
 * corners), the outernode data value belongs to only one data array so that 
 * there are no redundant data values.  Specifically, when DIM > 1, outernode 
 * data boxes are "trimmed" so that each node index that lives on more than one 
 * face on the box boundary will be associated with the face of the largest 
 * coordinate direction and only that face.  Within each array, data is stored 
 * in (i,...,k,d) order, where i,...,k indicates a spatial index and the d 
 * indicates the component depth at that location.  Memory allocation is
 * in column-major ordering (e.g., Fortran style) so that the leftmost
 * index runs fastest in memory.
 *
 * To illustrate the outernode data layout, in particular the "box trimming"
 * that prevents redundant data values, we describe the data for a
 * three-dimensional outernode data object instantiated over a box 
 * [l0:u0,l1:u1,l2:u2] in the standard SAMRAI cell-centered AMR index space.
 *
 * \verbatim

      Here face normal directions 0, 1, and 2 can be thought of as X, Y, Z
      respectively, and d is the data depth.

      face normal 0:
          lower    [ l0 : l0     , l1+1 : u1   , l2+1 : u2 , d ]
          upper    [ u0+1 : u0+1 , l1+1 : u1   , l2+1 : u2 , d ]
          Note: Boxes are trimmed at edges intersecting faces with
                normal directions 1 and 2 so that node indices shared 
                with those faces appear in data arrays associated with
                higher dimension faces.
          
      face normal 1:
          lower    [ l0 : u0+1   , l1 : l1     , l2+1 : u2 , d ]
          upper    [ l0 : u0+1   , u1+1 : u1+1 , l2+1 : u2 , d ]
          Note: Boxes are trimmed at edges intersecting faces with
                normal direction 2 so that node indices shared 
                with those faces appear in data arrays associated with
                higher dimension faces.
 
      face normal 2:
          lower    [ l0 : u0+1   , l1 : u1+1   , l2 : l2     , d ]
          upper    [ l0 : u0+1   , l1 : u1+1   , u2+1 : u2+1 , d ]
          Note: Boxes are not trimmed. 
  
 * \endverbatim
 * Other spatial dimensions are represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::OuternodeDataFactory
 * @see pdat::OuternodeGeometry
 * @see pdat::NodeIterator
 * @see pdat::NodeIndex
 */

template <int DIM, class TYPE>
class OuternodeData : public hier::PatchData<DIM>
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent outernode-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of
    * memory needed for TYPE is sizeof(TYPE).
    * If this is not the case, then a specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outernode data object will be created.
    *            Note: the ghost cell width is assumed to be zero.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    */
   static size_t getSizeOfData(const hier::Box<DIM>& box,
                               int depth);

   /*!
    * @brief Constructor for an outernode data object.
    *
    * Note: Outernode data always has ghost cell width of zero.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            outernode data object will be created.
    * @param depth gives the number of data values for each
    *              spatial location in the array.
    * @param pool memory arena.  If not given, then the
    *             standard arena is used.
    */
   OuternodeData(const hier::Box<DIM>& box,
                 int depth,
                 tbox::Pointer<tbox::Arena> pool = 
                    tbox::Pointer<tbox::Arena>(NULL));

   /*!
    * @brief Virtual destructor for a outernode data object.
    */
   virtual ~OuternodeData<DIM,TYPE>();

   /*!
    * @brief Return the depth (e.g., the number of components at each spatial
    * location) of the array.
    */
   int getDepth() const;

   /*!
    * @brief Returns true if outernode data exists for the given
    * face normal direction; false otherwise.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    */
   bool dataExists(int face_normal) const;

   /*!
    * @brief Return the box of valid node indices for
    *        outernode data.  Note: the returned box
    *        will reside in the @em node index space.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    */
   hier::Box<DIM> getDataBox(int face_normal,
                             int side);

   /*!
    * @brief Get a pointer to the beginning of a particular face normal, 
    * side, and depth component of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   TYPE* getPointer(int face_normal,
                    int side,
                    int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular face 
    * normal, side, and depth component of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   const TYPE* getPointer(int face_normal, 
                          int side,
                          int depth = 0) const;

   /*!
    * @brief Return a reference to data entry corresponding
    * to a given node index and depth.
    *
    * @param i const reference to NodeIndex, @em MUST be
    *          an index on the outernode of the box.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   TYPE& operator()(const NodeIndex<DIM>& i, 
                    int depth = 0);

   /*!
    * @brief Return a const reference to data entry corresponding
    * to a given node index and depth.
    *
    * @param i const reference to NodeIndex, @em MUST be
    *          an index on the outernode of the box.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    */
   const TYPE& operator()(const NodeIndex<DIM>& i, 
                          int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for
    * face normal, and side index of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    */
   ArrayData<DIM,TYPE>& getArrayData(int face_normal,
                                     int side);
 
   /*!
    * @brief Return a const reference to the array data object for
    * face normal, and side index of the outernode centered array.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outeredge
    *             data array
    */
   const ArrayData<DIM,TYPE>& getArrayData(int face_normal,
                                           int side) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * either NodeData or OuternodeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results. 
    */
   virtual void copy(const hier::PatchData<DIM>& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * either NodeData or OuternodeData of the same DIM and TYPE.  If not,
    * then an unrecoverable error results.
    */
   virtual void copy2(hier::PatchData<DIM>& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap. 
    *
    * Currently, source data must be either NodeData or OuternodeData 
    * of the same DIM and TYPE and the overlap must be an NodeOverlap
    * of the same DIM.  If not, then an unrecoverable error  
    * results. 
    */
   virtual void copy(const hier::PatchData<DIM>& src,
                     const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination 
    * patch data object on the given overlap.
    *
    * Currently, destination data must be either NodeData or OuternodeData  
    * of the same DIM and TYPE and the overlap must be an NodeOverlap
    * of the same DIM.  If not, then an unrecoverable error  
    * results. 
    */
   virtual void copy2(
      hier::PatchData<DIM>& dst,
      const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Fast copy (i.e., assumes node and outernode data objects are
    * defined over the same box) from the given node source data object to 
    * this destination outernode data object at the specified depths.
    */
   void copyDepth(int dst_depth,
                  const NodeData<DIM,TYPE>& src,
                  int src_depth);

   /*!
    * @brief Fast copy (i.e., assumes node and outernode data objects are
    * defined over the same box) to the given node destination data object 
    * from this source outernode data object at the specified depths.
    */
   void copyDepth2(int dst_depth,
                   NodeData<DIM,TYPE>& dst,
                   int src_depth) const;

   /*!
    * @brief Add data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be OuternodeData of the same DIM and 
    * TYPE and the overlap must be an EdgeOverlap of the same DIM.  
    * If not, then an unrecoverable error results.
    */
   virtual void sum(const hier::PatchData<DIM>& src,
                    const hier::BoxOverlap<DIM>& overlap);

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
    * @brief Pack data in this patch data object lying in the specified 
    * box overlap region into the stream.  The overlap must be an
    * NodeOverlap of the same DIM.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const hier::BoxOverlap<DIM>& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over 
    * the specified box overlap region.  The overlap must be an
    * NodeOverlap of the same DIM.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const hier::BoxOverlap<DIM>& overlap);

   /*!
    * @brief Unpack data from stream and add into this patch data object 
    * over the specified box overlap region.  The overlap must be an
    * NodeOverlap of the same DIM.
    */
   virtual void unpackStreamAndSum(tbox::AbstractStream& stream,
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
    * @brief Print all outernode data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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
    * @brief Print all outernode data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
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
    * @brief Print all outernode centered data values for specified 
    * face_normal and side residing in the specified box.
    * If the depth of the data is greater than one, all depths are printed.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void printAxisSide(int face_normal,
                      int side,
                      const hier::Box<DIM>& box,
                      std::ostream& os = tbox::plog,
                      int prec = 12) const;

   /*!
    * @brief Print all outernode centered data values for specified 
    * face_normal, side, and depth residing in the specified box.
    *
    * @param face_normal  integer face normal direction for data,
    *              must satisfy 0 <= face_normal < DIM
    * @param side integer lower (0) or upper (1) side of outernode
    *             data array
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to node index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void printAxisSide(int face_normal,
                      int side,
                      const hier::Box<DIM>& box,
                      int depth,
                      std::ostream& os = tbox::plog,
                      int prec = 12) const;

   /*!
    * @brief Check that class version and restart file version are equal.  
    * If so, read data members from the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void getSpecializedFromDatabase(
        tbox::Pointer<tbox::Database> database);

   /*!
    * @brief Write out the class version number and other data members to
    * the database.
    *
    * Assertions: database must be a non-null pointer.
    */
   virtual void putSpecializedToDatabase(
        tbox::Pointer<tbox::Database> database); 

private:
   OuternodeData<DIM,TYPE>(const OuternodeData<DIM,TYPE>&); // not implemented
   void operator=(const OuternodeData<DIM,TYPE>&);	  // not implemented

   //@
   //! @name Internal implementations of data copy operations.
   void copyFromNode( const NodeData<DIM,TYPE> &src );
   void copyFromNode( const NodeData<DIM,TYPE> &src,
		      const NodeOverlap<DIM> &overlap  );
   void copyToNode( NodeData<DIM,TYPE> &dst ) const;
   void copyToNode( NodeData<DIM,TYPE> &dst,
		    const NodeOverlap<DIM> &overlap  ) const;
   void copyFromOuternode( const OuternodeData<DIM,TYPE> &src );
   void copyFromOuternode( const OuternodeData<DIM,TYPE> &src,
			   const NodeOverlap<DIM> &overlap );
   void copyToOuternode( OuternodeData<DIM,TYPE> &dst ) const;
   void copyToOuternode( OuternodeData<DIM,TYPE> &dst,
			 const NodeOverlap<DIM> &overlap ) const;
   //@
		 

   int d_depth;
   ArrayData<DIM,TYPE> d_data[DIM][2];

};

}
}

#ifndef DEBUG_NO_INLINE
#include "OuternodeData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuternodeData.C"
#endif
