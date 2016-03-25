//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/CoarsenClasses.h $
// Package:	SAMRAI data transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Simple structure for managing coarsening data in equivalence classes.
//
 
#ifndef included_xfer_CoarsenClasses
#define included_xfer_CoarsenClasses

#include "SAMRAI_config.h"

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#include "tbox/Array.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#include "CoarsenOperator.h"

namespace SAMRAI {
    namespace xfer {

#ifndef NULL
#define NULL (0)
#endif

/**
 * Class CoarsenClasses<DIM> is used by the CoarsenSchedule<DIM> and 
 * CoarsenAlgorithm<DIM> classes to manage a set of data coarsening items 
 * that describe coarsening of patch data between two levels in an AMR hierarchy.  
 * Specifically, this class organizes these items into equivalence clases.
 * Two items are equivalent if their source and destination data represent 
 * the same data types and have the same ghost cell widths, respectively.
 */

template<int DIM> class CoarsenClasses : public tbox::DescribedClass
{
public:
   /**
    * The <TT>CoarsenClasses<DIM>::Data</TT> data structure contains the items that 
    * describe data coarsening between two levels on an AMR hierarchy.
    */
   struct Data {
      int                                         d_dst;
      int                                         d_src;
      bool                                        d_fine_bdry_reps_var;
      hier::IntVector<DIM>                        d_gcw_to_coarsen;
      tbox::Pointer< xfer::CoarsenOperator<DIM> > d_opcoarsen;
      int                                         d_tag;
   };

   /**
    * The constructor sets boolean for filling coarse data and creates 
    * an empty array of coarsen classes.
    */
   CoarsenClasses(bool fill_coarse_data);

   /**
    * The virtual destructor destroys the coarsen data items owned
    * by this object (and the associated CoarsenAlgorithm<DIM> object).
    */
   virtual ~CoarsenClasses<DIM>();

   /**
    * Return number of equivalence classes maintained by this object
    * (i.e., the number of lists of coarsen data items).
    */
   int getNumberOfEquivalenceClasses() const;

   /**
    * Return number of coarsen data items in the equivalence classes
    * represented by the given integer identifier.
    */
   int getNumberOfItemsInEquivalenceClass(int equiv_class_id) const;

   /**
    * Return const reference to representative element of equivalence class 
    * with the given integer identifier.  When assertion checking is active, 
    * the id will be checked for validity.
    */
   const typename CoarsenClasses<DIM>::Data&
      getClassRepresentative(int equiv_class_id) const;
 
   /**
    * Return an iterator for the coarsen data list corresponding to the
    * equivalence class with the given integer identifier.  The number of
    * equivalence classes can be determined via the
    * getNumberOfEquivalenceClasses() member function.  Valid integer
    * arguments are 0,...,getNumberOfEquivalenceClasses()-1.  When assertion 
    * checking is active, the id will be checked for validity.
    *
    * Note that the list should not be modified through this iterator.
    * When assertion checking is active, the id will be checked for validity.
    */
#ifdef LACKS_NAMESPACE_IN_DECLARE
   typename tbox::List<Data>::Iterator getIterator(int equiv_class_id);
#else
   typename tbox::List<typename CoarsenClasses<DIM>::Data>::Iterator
      getIterator(int equiv_class_id);
#endif

   /**
    * Insert a data item into the coarsen data list for the proper 
    * equivalence class.  Items are inserted in order of operator
    * priority so that communication algorithms apply coarsen operators
    * with the lowest numerical priority before those with higher numerical
    * priority.
    *
    * If a null patch descriptor argument is passed (or ommitted), the descriptor
    * associated with the variable database Singleton object will be used to
    * determine the equivalence class.
    */
   void insertEquivalenceClassItem(
      const typename CoarsenClasses<DIM>::Data& data, 
      tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                        (hier::PatchDescriptor<DIM>*)NULL);

   /**
    * Check coarsen data item so that the source and destination
    * data entries have an appropriate number of ghost cells to
    * support the stencil width of the operator passed to the 
    * CoarsenAlgorithm::registerCoarsen() routine.  Also, if 
    * the boolean argument passed to CoarsenAlgorithm constructor 
    * was true, indicating that coarse level data must be filled before 
    * coarsening, this routine checks that copying between source and 
    * destination data types is a valid operation.  If both checks are 
    * successful, return true; else return false.  A descriptive error 
    * message will be reported when a problem appears and the program halts.
    *
    * If a null patch descriptor argument is passed (or omitted), the 
    * descriptor ownded by the variable database Singleton object is used.
    */
   bool checkCoarsenItem(const typename CoarsenClasses<DIM>::Data& data_item,
                         tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                         (hier::PatchDescriptor<DIM>*)NULL) const;

   /**
    * Compare coarsen data items in this coarsen classes object against those
    * in the argument coarsen classes object.  Return true if they all match
    * with regard to the patch data types, patch data ghost cell widths,
    * operator stencils, etc. that they refer to and return false otherwise.
    *
    * If a null patch descriptor argument is passed (or ommitted), the descriptor
    * associated with the variable database Singleton object will be used.
    */
   bool checkConsistency(tbox::Pointer< CoarsenClasses<DIM> > test_classes,
                         tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                         (hier::PatchDescriptor<DIM>*)NULL) const;

   /**
    * Print all equivalence class data to the specified output stream.
    */
   virtual void printClassData(std::ostream& stream) const;

   /**
    * Print single equivalence class item to the specified output stream.
    */
   void printCoarsenItem(std::ostream& stream,
                        const typename CoarsenClasses<DIM>::Data& data) const;

private:
   CoarsenClasses(const CoarsenClasses<DIM>&);	// not implemented
   void operator=(const CoarsenClasses<DIM>&);         // not implemented

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
      const typename CoarsenClasses<DIM>::Data& data,
      tbox::Pointer< hier::PatchDescriptor<DIM> > descriptor =
                   (hier::PatchDescriptor<DIM>*)NULL) const;

   /*
    * Boolean set in constructor to indicate whether coarse data will be
    * filled before coarsening operations are performed.
    */
   bool d_fill_coarse_data;

#ifdef LACKS_NAMESPACE_IN_DECLARE
   tbox::Array< tbox::List<Data> > d_coarsen_equivalence_classes;
#else
   tbox::Array< tbox::List<typename CoarsenClasses<DIM>::Data> > 
      d_coarsen_equivalence_classes;
#endif

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarsenClasses.C"
#endif
