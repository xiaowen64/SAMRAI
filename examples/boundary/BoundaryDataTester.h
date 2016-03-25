//
// File:        BoundaryDataTester.h
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Class to test usage of boundary utilities
//

#ifndef included_BoundaryDataTester
#define included_BoundaryDataTester

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_appu_BoundaryUtilityStrategy
#include "BoundaryUtilityStrategy.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_geom_CartesianGridGeometry
#include "CartesianGridGeometry.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_hier_VariableContext
#include "VariableContext.h"
#endif

using namespace SAMRAI;

class BoundaryDataTester : 
   public xfer::RefinePatchStrategy<NDIM>,
   public appu::BoundaryUtilityStrategy
{
public:
  /**
   * The constructor reads variable data from input database.
   */
   BoundaryDataTester(const string& object_name,
                      tbox::Pointer<tbox::Database> input_db,
                      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geom);

   /**
    * Virtual destructor for BoundaryDataTester.
    */
   virtual ~BoundaryDataTester();

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class RefinePatchStrategy.  It sets the boundary 
    * conditions for the variables. 
    */
   void setPhysicalBoundaryConditions(hier::Patch<NDIM>& patch,
                                      const double fill_time,
                                      const hier::IntVector<NDIM>& ghost_width_to_fill);

   /**
    * The next three functions are dummy implementations of the pure
    * virtual functions declared in the RefinePatchStrategy base class.
    * They are not needed for this example since we only have one level
    * in the hierarchy.
    */
   hier::IntVector<NDIM> getRefineOpStencilWidth() const { return hier::IntVector<NDIM>(0); }

   void preprocessRefine(hier::Patch<NDIM>& fine,
                         const hier::Patch<NDIM>& coarse,
                         const hier::Box<NDIM>& fine_box,
                         const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   }

   void postprocessRefine(hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio)
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   }

   /**
    * This routine is a concrete implementation of a virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * face or edge boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face or edge to which the boundary condition applies.
    */
   void readDirichletBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                       string& db_name,
                                       int bdry_location_index);

   /**
    * This routine is a concrete implementation of a virtual function
    * in the base class BoundaryUtilityStrategy.  It reads NEUMANN
    * face or edge boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face or edge to which the boundary condition applies.
    */
   void readNeumannBoundaryDataEntry(tbox::Pointer<tbox::Database> db,
                                     string& db_name,
                                     int bdry_location_index);

   /**
    * Set data on patch interiors on given level in hierarchy.
    */
   void initializeDataOnPatchInteriors(tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                       int level_number);

   /**
    * Run boundary tests for given level in hierarchy.
    */
   void runBoundaryTest(tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                        int level_number);

   /**
    * Print all class data members to given output stream.
    */
   void printClassData(ostream &os) const;

private:
   /*
    * The object name is used for error/warning reporting.
    */
   string d_object_name;

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

   /*
    * Arrays of information read from input file describing test variables
    */
   tbox::Array<string>               d_variable_name;
   tbox::Array<int>                  d_variable_depth;
   tbox::Array<hier::IntVector<NDIM> >            d_variable_num_ghosts;
   tbox::Array< tbox::Array<double> > d_variable_interior_values;

   /*
    * Items used to manage variables and data in test program.
    */
   tbox::Array< tbox::Pointer<hier::Variable<NDIM> > >  d_variables;
   tbox::Pointer<hier::VariableContext> d_variable_context;
   hier::ComponentSelector d_patch_data_components;

   /*
    * Arrays of information read from input file for boundary conditions
    */ 
   tbox::Array<int> d_master_bdry_edge_conds;
   tbox::Array<int> d_scalar_bdry_edge_conds; 
   tbox::Array<int> d_vector_bdry_edge_conds; 

   tbox::Array<int> d_master_bdry_node_conds;
   tbox::Array<int> d_scalar_bdry_node_conds; 
   tbox::Array<int> d_vector_bdry_node_conds; 

#if (NDIM == 3)
   tbox::Array<int> d_master_bdry_face_conds;
   tbox::Array<int> d_scalar_bdry_face_conds;
   tbox::Array<int> d_vector_bdry_face_conds;
#endif

#if (NDIM == 2)
   tbox::Array<int> d_node_bdry_edge;
#endif
#if (NDIM == 3)
   tbox::Array<int> d_edge_bdry_face;
   tbox::Array<int> d_node_bdry_face;
#endif

   tbox::Array< tbox::Array<double> > d_variable_bc_values;

   /*
    * Private functions to perform tasks for boundary testing.
    */
   void readVariableInputAndMakeVariables(tbox::Pointer<tbox::Database> db);
   void readBoundaryDataInput(tbox::Pointer<tbox::Database> db);
   void readBoundaryDataStateEntry(tbox::Pointer<tbox::Database> db,
                                   string& db_name,
                                   int bdry_location_index); 
   void setBoundaryDataDefaults();
   void postprocessBoundaryInput();
   void checkBoundaryData(int btype,
                          const hier::Patch<NDIM>& patch,
                          const hier::IntVector<NDIM>& ghost_width_to_check) const;

};

#endif
