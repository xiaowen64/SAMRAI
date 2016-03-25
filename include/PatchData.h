//
// File:	PatchData.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Abstract base class for patch data objects
//

#ifndef included_hier_PatchData
#define included_hier_PatchData

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_AbstractStream
#include "tbox/AbstractStream.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxOverlap
#include "BoxOverlap.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif


namespace SAMRAI {
    namespace hier {

/**
 * Class PatchData<DIM> is a pure virtual base class for the data storage
 * defined over a box.  Patch data objects are generally contained within
 * a patch.  Patch data defines the abstract virtual functions for data
 * management and communication that must be supplied by the subclasses.
 * Subclasses implement the virtual functions as appropriate for the derived
 * type; for example, cell centered objects will copy data differently than
 * face centered objects.
 *
 * Patch data objects are created by patch data factories and associated
 * subclasses.  This separation into abstract factory and concrete
 * implementation subclasses facilitates the creation of new patch data
 * subtypes.  See the Design Patterns book for more details about the
 * Abstract Factory pattern.
 *
 * The copy and pack/unpack functions in the patch data object take
 * box overlap descriptions that describe the index space over which
 * data is to be copied or packed/unpacked.  Box overlaps are computed
 * by the box geometry classes, which are accessed through the patch
 * data factories.  Box geometry classes are created by the patch data
 * factories instead of the patch data objects because patch data objects
 * are distributed across memory with patches and therefore may not exist
 * on a particular processor.  Patch data factories are guaranteed to
 * exist on all processors independent of the patch-to-processor mapping.
 *
 * @see hier::BoxOverlap
 * @see hier::BoxGeometry
 * @see hier::Patch
 * @see hier::PatchDataFactory
 * @see hier::PatchDescriptor
 */

template<int DIM> class PatchData : public tbox::DescribedClass
{
public:
   /**
    * The constructor for a patch data object.  Patch data objects will
    * manage the interior box over which they are defined and the associated
    * ghost cell width.
    */
   PatchData(const Box<DIM>& domain, const IntVector<DIM>& ghosts);

   /**
    * The virtual destructor for a patch data object.
    */
   virtual ~PatchData<DIM>();

   /**
    * Return the box over which this patch data object is defined.  All
    * objects in the same patch are defined over the same box, although
    * the patch data objects may interpret how to allocate storage for
    * that box in different ways.
    */
   const Box<DIM>& getBox() const;

   /**
    * Return the ghost cell box.  The ghost cell box is defined to be
    * the interior box grown by the ghost cell width.
    */
   const Box<DIM>& getGhostBox() const;

   /**
    * Get the ghost cell width associated with this patch data object.
    */
   const IntVector<DIM>& getGhostCellWidth() const;

   /**
    * Set the simulation time stamp for the patch data type.  The simulation
    * time is initialized to zero when the patch data type is created.
    */
    void setTime(const double timestamp);

   /**
    * Get the simulation time stamp for the patch data type.
    */
   double getTime() const;

   /**
    * A fast copy between the source and destination.  Data is copied from
    * the source into the destination where there is overlap in the underlying
    * index space.  The copy is performed on the interior plus the ghost cell
    * width (for both the source and destination).  If this copy does not
    * understand how to copy data from the argument, then copy2() is called
    * on the source object.
    */
   virtual void copy(const PatchData<DIM>& src) = 0;

   /**
    * A fast copy between the source and destination.  Data is copied from
    * the source into the destination where there is overlap in the underlying
    * index space.  The copy is performed on the interior plus the ghost cell
    * width (for both the source and destination).  If this copy does not
    * understand how to copy data from the destination, then it may throw
    * an exception (aka dump core in a failed assertion).
    */
   virtual void copy2(PatchData<DIM>& dst) const = 0;

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects.  If this member function
    * cannot complete the copy from source (e.g., if it doesn't understand
    * the type of source), then copy2() is called on the source object.
    */
   virtual void copy(const PatchData<DIM>& src,
                     const BoxOverlap<DIM>& overlap) = 0;

   /**
    * Copy data from the source into the destination using the designated
    * overlap descriptor.  The overlap description will have been computed
    * using the appropriate box geometry objects If this member function
    * cannot complete the copy from the destination, then it may throw an
    * exception (aka dump core in a failed assertion).
    */
   virtual void copy2(PatchData<DIM>& dst,
                      const BoxOverlap<DIM>& overlap) const = 0;

   /**
    * Determines whether the patch data subclass can estimate the necessary
    * stream size using only index space information.  The return value will
    * most likely be true for data types that are fixed size (such as doubles)
    * but will be false for complex data types that allocate their own storage
    * (such as lists of particles).  This routine is used to estimate whether
    * a processor can estimate space for incoming messages or whether it needs
    * to receive a message size from the sending processor.
    */
   virtual bool canEstimateStreamSizeFromBox() const = 0;

   /**
    * Calculate the number of bytes needed to stream the data lying
    * in the specified box domain.  This estimate must be an upper
    * bound on the size of the data in the actual message stream.
    * The upper bound should be close, however, since buffer space
    * will be allocated according to these values, and excess buffer
    * space will waste memory resources.
    */
   virtual int getDataStreamSize(const BoxOverlap<DIM>& overlap) const = 0;

   /**
    * Pack data lying on the specified index set into the output stream.
    * See the abstract stream virtual base class for more information about
    * the packing operators defined for streams.
    */
   virtual void packStream(tbox::AbstractStream& stream,
                           const BoxOverlap<DIM>& overlap) const = 0;

   /**
    * Unpack data from the message stream into the specified index set.
    * See the abstract stream virtual base class for more information about
    * the packing operators defined for streams.
    */
   virtual void unpackStream(tbox::AbstractStream& stream,
                             const BoxOverlap<DIM>& overlap) = 0;

   /**
    * Checks that class version and restart file version are equal.  If so,
    * reads in the data members common to all patch data types from database.
    * This method then calls the getSpecializedFromDatabase() method
    * to retrieve the data special to the concrete patch data type.
    */
   virtual void getFromDatabase(tbox::Pointer<tbox::Database> database);

   /**
    * Writes out the class version number to the database.  Then, 
    * writes the data members common to all patch data types from database.
    * After the common data is written to the database, the 
    * putSpecializedToDatabase() method is invoked.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> database);

   /**
    * This pure abstract method is used by concrete patch data subclasses 
    * to retrieve from the database data special to the concrete class.
    */
   virtual void getSpecializedFromDatabase(
           tbox::Pointer<tbox::Database> database) = 0;

   /**
    * This pure abstract method is used by concrete patch data subclasses 
    * to put to the database data special to the concrete class.
    */
   virtual void putSpecializedToDatabase(
           tbox::Pointer<tbox::Database> database) = 0;

protected:
   /**
    * This protected method is used by concrete patch data subclasses
    * to set the ghost box over which the patch data will be allocated.
    * Note that this allows the ghost box to be inconsistant with its
    * standard interpretation as the patch domain box grown by the ghost
    * cell width (as set in the constructor).  
    *
    * This function is included to treat some special cases for concrete
    * patch data types and should be used with caution.
    */
   void setGhostBox(const Box<DIM>& ghost_box);

private:
   PatchData(const PatchData<DIM>&);	// not implemented
   void operator=(const PatchData<DIM>&);	// not implemented

   Box<DIM> d_box;				// interior box description
   Box<DIM> d_ghost_box;			// interior box plus ghosts
   IntVector<DIM> d_ghosts;			// ghost cell width
   double d_timestamp;				// timestamp for the data

};

}
}
#ifndef DEBUG_NO_INLINE
#include "PatchData.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "PatchData.C"
#endif
