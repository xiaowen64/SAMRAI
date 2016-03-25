//
// File:	Variable.h
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 301 $
// Modified:	$Date: 2005-04-25 10:32:08 -0700 (Mon, 25 Apr 2005) $
// Description:	Base class for application-level variables
//

#ifndef included_hier_Variable
#define included_hier_Variable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_PatchDataFactory
#include "PatchDataFactory.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif

namespace SAMRAI {
    namespace hier {

/**
 * Class Variable<DIM> provides a description of a variable quantity in 
 * an AMR application.  A variable has the following attributes: (1) a name, 
 * (2) a unique instance number, and (3) a way to create the proper type of 
 * storage on the patches in an AMR hierarchy.  
 *
 * Variable is a base class for all variables; subclasses for each data type 
 * (e.g., cell-centered variables or node-centered variables) will define 
 * member functions to create the appropriate storage management classes.  
 * Note that a variable object is distinct from the storage for the quantity 
 * to which the variable refers.
 *
 * The user defines a set of variables and the role they play in the solution 
 * algorithm when describing an AMR application.  Each variable type knows
 * how to create instances of the associated data objects on patches in an
 * AMR patch hierarchy via the getPatchDataFactory() member function.  For 
 * example, we might have a concrete  myVariable3<double> object that is 
 * derived from Variable3.  The myVariable3<double> object returns a 
 * myVariableFactory3<double> object from its member function 
 * getPatchDataFactory(), and this factory is used to create the 
 * myVariable3<double> objects on the patches.  Thus, a solution algorithm 
 * may be implemented to manage storage of variable quantities through the 
 * abstract interface without knowing the specific concrete variable types
 * in use.  This approach gives maximum flexibility when defining new user 
 * data types and variables and solution algorithms.
 *
 * Each variable is assigned an ``instance number,'' a unique integer
 * identifier numbered from 0 through MAX-1, where MAX is the value
 * returned by static function getCurrentMaximumInstanceNumber().
 * These identifiers can be used to rapidly look up variable instances
 * in a table.  For example, the instance number can be used to map
 * variables to algorithmic data within an algorithm as alluded to above..
 *
 * Variable types for which data exists on the borders of patches (such as node,
 * side, etc.), the data will thus live on the interface between coarse and fine 
 * patch levels.  Thus, it must be specified whether coarse or fine data values 
 * take precedence on the interface between levels.  This information is provided 
 * by the fineBoundaryRepresentsVariable() function.  Each concrete variable subclass
 * defines the behavior of this function.
 * 
 * @see hier::PatchDataFactory
 */

template<int DIM> class Variable  : public tbox::DescribedClass
{
public:
   /**
    * Return the current maximum instance number over all variable objects.
    * The instance identifier returned from each variable objhect is 
    * guaranteed to be between 0 and this number minus one.  Note that this 
    * number changes as new variable instances are created.
    */
   static int getCurrentMaximumInstanceNumber();

   /**
    * Create a variable object with the specified name and patch data
    * factory.  On creation, each variable is assigned a unique instance
    * identifier.
    */
   Variable(const string &name,
                  const tbox::Pointer< PatchDataFactory<DIM> > factory);

   /**
    * Virtual destructor for variable objects.
    */
   virtual ~Variable<DIM>();

   /**
    * Return the instance identifier for this particular variable object.
    * The instance identifiers are unique integers numbered starting from zero.
    */
   int getInstanceIdentifier() const;

   /**
    * Return the name assigned to this variable.
    */
   const string& getName() const;

   /**
    * Return true if the fine data values represent the variable quantity
    * on coarse-fine interfaces if variable data lives on patch borders; 
    * false otherwise.  The boolean return value is supplied by the concrete 
    * variable subclass. 
    */
   virtual bool fineBoundaryRepresentsVariable() const = 0;

   /**
    * Return true if the variable data lives on patch borders; false otherwise.  
    * The boolean return value is supplied by the concrete variable subclass.
    */
   virtual bool dataLivesOnPatchBorder() const = 0;

   /**
    * Set the patch data factory object.  Normally, the factory is set in
    * the constructor, but this member function enables the factory to be
    * changed later in the lifetime of the variable.
    */
   void setPatchDataFactory(tbox::Pointer< PatchDataFactory<DIM> > factory);

   /**
    * Return a non-const pointer to a patch data factory that will be used
    * to instantiate instances of this variable on the patches.  The factory
    * returned will have been set by the variable subclasses.
    */
   tbox::Pointer< PatchDataFactory<DIM> > getPatchDataFactory() const;

private:
   Variable(const Variable<DIM>&);	// not implemented
   void operator=(const Variable<DIM>&);	// not implemented

   string d_name;
   int d_instance;
   tbox::Pointer< PatchDataFactory<DIM> > d_factory;

   static int s_instance_counter;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "Variable.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "Variable.C"
#endif
