//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/RefineClasses.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Simple structure for managing refinement data in equivalence classes.
//
 
#ifndef included_xfer_RefineClasses
#define included_xfer_RefineClasses

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#include "RefineOperator.h"
#include "TimeInterpolateOperator.h"

namespace SAMRAI {
    namespace xfer {

/**
 * Class RefineClasses<DIM> is used by the RefineSchedule<DIM> and 
 * RefineAlgorithm<DIM> classes to manage a set of refinement data items 
 * that describe interpatch communication of patch data on an AMR hierarchy.  
 * Specifically, this class organizes these items into equivalence clases.
 * Two items are equivalent if their source and destination data represent 
 * the same data types and have the same ghost cell widths, respectively.
 */

template<int DIM> class RefineClasses : public tbox::DescribedClass
{
public:
   /**
    * The <TT>RefineClasses<DIM>::Data</TT> data structure contains the items that 
    * describe communication between patch data components on an AMR hierarchy.
    */
   struct Data {
      int                                           d_dst;
      int                                           d_src;
      int                                           d_src_told;
      int                                           d_src_tnew;
      int                                           d_scratch;
      bool                                          d_fine_bdry_reps_var;
      bool                                          d_time_interpolate;
      tbox::Pointer< RefineOperator<DIM> >          d_oprefine;
      tbox::Pointer< TimeInterpolateOperator<DIM> > d_optime;
      int                                           d_tag;
   };

   /**
    * The constructor creates an empty array of refine classes.
    */
   RefineClasses();

   /**
    * The virtual destructor destroys the refinement data items owned 
    * by this object (and the associated RefineAlgorithm<DIM> object).
    */
   virtual ~RefineClasses<DIM>();

   /**
    * Return number of equivalence classes maintained by this object
    * (i.e., the number of lists of refinement data items).
    */
   int getNumberOfEquivalenceClasses() const;

   /**
    * Return const reference to representative element of equivalence class 
    * with the given integer identifier.  When assertion checking is active, 
    * the id will be checked for validity.
    */
   const typename RefineClasses<DIM>::Data&
      getClassRepresentative(int equiv_class_id) const;
 
   /**
    * Return an iterator for the refine data list corresponding to the
    * equivalence class with the given integer identifier.  The number of
    * equivalence classes can be determined via the
    * getNumberOfEquivalenceClasses() member function.  Valid integer
    * arguments are 0,...,getNumberOfEquivalenceClasses()-1. When assertion
    * checking is active, the id will be checked for validity.
    *
    * Note that the list should not be modified through this iterator.
    * When assertion checking is active, the id will be checked for validity.
    */
#ifdef LACKS_NAMESPACE_IN_DECLARE
   typename tbox::List<Data>::Iterator getIterator(int equiv_class_id);
#else
   typename tbox::List<typename RefineClasses<DIM>::Data>::Iterator
      getIterator(int equiv_class_id);
#endif

   /**
    * Insert a data item into the refine data list for the proper 
    * equivalence class.  Items are inserted in order of operator
    * priority so that communication algorithms apply refinement operators
    * with the lowest numerical priority before those with higher numerical
    * priority.
    *
    * If a null patch descriptor argument is passed (or ommitted), the descriptor
    * associated with the variable database Singleton object will be used to 
    * determine the equivalence class.
    */
   void insertEquivalenceClassItem(
      const typename RefineClasses<DIM>::Data& data, 
      tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                        (hier::PatchDescriptor<DIM>*)NULL);

   /**
    * Check refine data item so that scratch data entry has at least as
    * many ghost cells as destination data entry and stencil width of
    * operator.  If so, return true; else return false.  A descriptive
    * error message will be reported when a problem appears and the 
    * program halts.
    *
    * If a null patch descriptor argument is passed (or ommitted), the descriptor 
    * associated with the variable database Singleton object will be used.
    */
   bool checkRefineItem(const typename RefineClasses<DIM>::Data& data_item,
                        tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                        (hier::PatchDescriptor<DIM>*)NULL) const;

   /**
    * Compare refine data items in this refine classes object against those
    * in the argument refine classes object.  Return true if they all match 
    * with regard to the patch data types, patch data ghost cell widths,
    * operator stencils, etc. that they refer to and return false otherwise.  
    *
    * If a null patch descriptor argument is passed (or ommitted), the descriptor
    * associated with the variable database Singleton object will be used.
    */
   bool checkConsistency(tbox::Pointer< RefineClasses<DIM> > test_classes,
                         tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                         (hier::PatchDescriptor<DIM>*)NULL) const;

   /**
    * Print all equivalence class data to the specified output stream.
    */
   virtual void printClassData(std::ostream& stream) const;

   /**
    * Print single equivalence class item to the specified output stream.
    */
   void printRefineItem(std::ostream& stream,
                        const typename RefineClasses<DIM>::Data& data) const;

private:
   RefineClasses(const RefineClasses<DIM>&);	// not implemented
   void operator=(const RefineClasses<DIM>&);         // not implemented

   /*
    * Function to compare two patch data items (with given descriptor indices)
    * for consistency.  Return true if consistent; false otherwise.
    */
   bool checkPatchDataItemConsistency(
      int item_id1,
      int item_id2,
      tbox::Pointer< hier::PatchDescriptor<DIM> > pd) const;

   /*
    * Function to compare destination and source pair in data item to previously
    * registered pairs, to determine if they are equivalent to an operation
    * already registered with this algorithm.  Equivalence is defined as
    * representing the same data type and having the same ghost cell widths.
    * If the operation is equivalence to an equivalence class that has
    * been previously found, an integer identifier for that equivalence
    * class is returned.  If not it returns -1, an invalid equivalence index.
    */
   int getEquivalenceClassIndex(
      const typename RefineClasses<DIM>::Data& data,
      tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                   (hier::PatchDescriptor<DIM>*)NULL) const;

#ifdef LACKS_NAMESPACE_IN_DECLARE
   tbox::Array< tbox::List<Data> > d_refine_equivalence_classes;
#else
   tbox::Array< tbox::List<typename RefineClasses<DIM>::Data> > 
      d_refine_equivalence_classes;
#endif

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "RefineClasses.C"
#endif
