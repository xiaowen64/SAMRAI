//
// File:	LocallyActiveDataFillBox.h
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 481 $
// Modified:	$Date: 2005-07-21 13:50:43 -0700 (Thu, 21 Jul 2005) $
// Description:	Routines for "smart" boxlist ops in locally-active data comm ops
//

#ifndef included_hier_LocallyActiveDataFillBox
#define included_hier_LocallyActiveDataFillBox

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_tbox_List
#include "tbox/List.h"
#endif
#ifndef included_tbox_PIO
#include "tbox/PIO.h"
#endif
#ifndef included_xfer_CoarsenClasses
#include "CoarsenClasses.h"
#endif
#ifndef included_xfer_RefineClasses
#include "RefineClasses.h"
#endif
#ifndef included_iostream
#include <iostream>
using namespace std;
#define included_iostream
#endif

namespace SAMRAI {
   namespace xfer {

/*!
 * Class LocallyActiveDataFillBox is a utility class that is used primarily
 * by the LocallyActiveDataFillBoxSet class for boxlist operations in communication 
 * schedules that operate on "locally-active" data; i.e., where each data item may 
 * live on a different set of patches.  Specifically, this class contains a box and 
 * an associated list of either CoarseClass or RefineClass items, but not both.  
 * Each constructor accepts such a list. Once an object is constructed, it can only 
 * be used to support either refine operations or coarsen operations.
 *
 * @see hier::Box
 * @see xfer::CoarsenClasses
 * @see xfer::RefineClasses
 */

template<int DIM>
class LocallyActiveDataFillBox
{
public:
   /*!
    * @brief Construct a locally-active fill box with the given box
    * and list of refine items.  Note that this object maintains pointers
    * to the given box and refine items only.  So, an object of this class
    * becomes unusable if those items are destroyed before this object.
    *
    * @param box       Input box.
    * @param var_data  Input list of refine class data pointers.
    */
   LocallyActiveDataFillBox(
      const hier::Box<DIM>& box,
      const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data);

   /*!
    * @brief Construct a locally-active fill box with the given box
    * and list of coarsen items.  Note that this object maintains pointers
    * to the given box and coarsen items only.  So, an object of this class
    * becomes unusable if those items are destroyed before this object.
    *
    * @param box       Input box.
    * @param var_data  Input list of coarsen class data pointers.
    */
   LocallyActiveDataFillBox(
      const hier::Box<DIM>& box,
      const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data);

   /*!
    * Copy constructor to create a new locally-active fill box and copy 
    * the information from the argument locally-active fill box.
    * 
    * @param fill_box  Constant reference to fill box to be copied.
    */
   LocallyActiveDataFillBox(const LocallyActiveDataFillBox<DIM>& fill_box);

   /*!
    * The destructor releases all box and locally-active data storage.
    */
   virtual ~LocallyActiveDataFillBox<DIM>();

   /*!
    * Return constant reference to box maintained by this object.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * Return constant reference to list of refine items maintained by this object.
    * 
    * Note that if this object was created using coarsen item data, an unrecoverable
    * assertion will result.
    */
   const tbox::List<const typename RefineClasses<DIM>::Data*>& 
      getActiveRefineVarData() const;

   /*!
    * Return constant reference to list of coarsen items maintained by this object.
    * 
    * Note that if this object was created using refine item data, an unrecoverable
    * assertion will result.
    */
   const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& 
      getActiveCoarsenVarData() const;

   /*!
    * Print all class member data for this locally-active fill box object
    * to specified output stream.
    */
   void printClassData(ostream& os) const;

   /*!
    * Clear all class member data for this locally-active fill box object.
    */
   void clearLocallyActiveFillBoxData(); 

   /*!
    * @brief Check given box and list of refine items for equality with those maintained
    * by this class object.  Any encountered inequality will be sent to the given output
    * stream.
    * 
    * @return   Boolean true if equal, false otherwise.
    * @param    box        Input box for comparison.
    * @param    var_data   Input list of refine items for comparison.
    * @param    os         Input ostream for reporting mismatches.
    * 
    * Note that if this object was created using coarsen item data, an unrecoverable
    * assertion will result.
    */
   bool checkData(const hier::Box<DIM>& box,
                  const tbox::List<const typename xfer::RefineClasses<DIM>::Data*>& var_data,
                  ostream& os) const;

   /*!
    * @brief Check given box and list of coarsen items for equality with those maintained
    * by this class object.  Any encountered inequality will be sent to the given output
    * stream.
    * 
    * @return   Boolean true if equal, false otherwise.
    * @param    box        Input box for comparison.
    * @param    var_data   Input list of coarsen items for comparison.
    * @param    os         Input ostream for reporting mismatches.
    *
    * Note that if this object was created using refine item data, an unrecoverable
    * assertion will result.
    */
   bool checkData(const hier::Box<DIM>& box,
                  const tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*>& var_data,
                  ostream& os) const;

private:
   const hier::Box<DIM>* d_box;

   tbox::List<const typename xfer::RefineClasses<DIM>::Data*> d_refine_var_data; 
   tbox::List<const typename xfer::CoarsenClasses<DIM>::Data*> d_coarsen_var_data; 

   bool d_refine_data;
};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataFillBox.C"
#endif
