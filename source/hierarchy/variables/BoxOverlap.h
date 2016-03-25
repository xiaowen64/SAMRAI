//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/variables/BoxOverlap.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Base class that describes intersections between AMR boxes
//

#ifndef included_hier_BoxOverlap
#define included_hier_BoxOverlap

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif


namespace SAMRAI {
    namespace hier {

/**
 * Class BoxOverlap<DIM> is an abstract base class used to represent the
 * intersections between two AMR boxes.  The exact form of the overlap
 * is determined by the particular geometry implemented by the subclass.
 * For example, the rules for intersecting face centered boxes is very
 * different from that for cell centered boxes.
 *
 * The BoxOverlap<DIM> class provides two functions.  First, it serves
 * as a base class that can answer the question whether an intersection
 * is empty, and is therefore useful for determining communication
 * dependencies.  Second, it is a storage location for the exact form of
 * the intersection of the data residing on two boxes, which can be quite
 * complicated (for example, for face centered boxes).  In the second case,
 * access to the intersection data is via narrowing the interface via
 * type-safe type casting and using the subclass member functions.
 *
 * @see hier::BoxGeometry
 */

template<int DIM> class BoxOverlap  : public tbox::DescribedClass
{
public:
   /**
    * The default constructor for BoxOverlap<DIM> does nothing interesting.
    */
   BoxOverlap();

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~BoxOverlap<DIM>();

   /**
    * Return true if intersecting boxes have no communication dependencies.
    * Note that two boxes may communicate even if they do not intersect in
    * the underlying AMR index space (e.g., if data values exist at the
    * outside portions of the cells).
    */
   virtual bool isOverlapEmpty() const = 0;

   /**
    * Return the offset between the destination and source index spaces.
    * The destination index space is the source index space shifted 
    * by this amount.
    */
   virtual const IntVector<DIM>& getSourceOffset() const = 0;

   /**
    * Output the box overlap.
    */   
   virtual void print(std::ostream& os) const; 


private:
   BoxOverlap(const BoxOverlap<DIM>&);	// not implemented
   void operator=(const BoxOverlap<DIM>&);	// not implemented

};

}
}
#ifndef DEBUG_NO_INLINE
#include "BoxOverlap.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "BoxOverlap.C"
#endif
