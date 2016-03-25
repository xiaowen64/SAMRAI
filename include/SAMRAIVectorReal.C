//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/vectors/SAMRAIVectorReal.C $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Vector class for data on SAMRAI hierarchy.
//

#ifndef included_solv_SAMRAIVectorReal_C
#define included_solv_SAMRAIVectorReal_C

#include "SAMRAIVectorReal.h"

#include <typeinfo>
#include <float.h>
#include <math.h>

#include "tbox/MathUtilities.h"


#include "HierarchyCellDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "PatchLevel.h"
#include "VariableDatabase.h"
#include "CellVariable.h"
#include "EdgeVariable.h"
#include "FaceVariable.h"
#include "NodeVariable.h"
#include "SideVariable.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

#ifdef DEBUG_NO_INLINE
#include "SAMRAIVectorReal.I"
#endif
namespace SAMRAI {
    namespace solv {

#define DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE (10)


/*
*************************************************************************
*                                                                       *
* Initialize the static operators and counters.                         *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE> int SAMRAIVectorReal<DIM,TYPE>::s_instance_counter = 0;

template<int DIM, class TYPE> tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> >
   SAMRAIVectorReal<DIM,TYPE>::s_cell_ops = NULL;
template<int DIM, class TYPE> tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> >
   SAMRAIVectorReal<DIM,TYPE>::s_edge_ops = NULL;
template<int DIM, class TYPE> tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> >
   SAMRAIVectorReal<DIM,TYPE>::s_face_ops = NULL;
template<int DIM, class TYPE> tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> >
   SAMRAIVectorReal<DIM,TYPE>::s_node_ops = NULL;
template<int DIM, class TYPE> tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> >
   SAMRAIVectorReal<DIM,TYPE>::s_side_ops = NULL;



/*
*************************************************************************
*                                                                       *
* The constructor for SAMRAIVectorReal<DIM> objects initializes        *
* vector structure.                                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
SAMRAIVectorReal<DIM,TYPE>::SAMRAIVectorReal(
   const std::string &name,
   tbox::Pointer< hier::PatchHierarchy<DIM> > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT(   (coarsest_level >= 0)
          && (finest_level >= coarsest_level)
          && (finest_level <= hierarchy->getFinestLevelNumber()) );
#endif

   SAMRAIVectorReal<DIM,TYPE>::s_instance_counter++;

   if (name.empty()) {
      d_vector_name = "SAMRAIVectorReal";
   } else {
      d_vector_name = name;
   }
   d_hierarchy = hierarchy;
   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level; 

   d_number_components = 0;

   // Set default output stream
   d_output_stream = &tbox::plog;
}


/*
*************************************************************************
*                                                                       *
* Destructor for SAMRAIVectorReal<DIM>.                                *
* Component data storage is not deallocated here.                       *
*                                                                       *
*************************************************************************
*/
template<int DIM, class TYPE>
SAMRAIVectorReal<DIM,TYPE>::~SAMRAIVectorReal()
{

   SAMRAIVectorReal<DIM,TYPE>::s_instance_counter--;

   d_number_components = 0;

   d_component_variable.resizeArray(0); 
   d_component_data_id.resizeArray(0);
   d_component_operations.resizeArray(0);
   d_control_volume_data_id.resizeArray(0);

   d_variableid_2_vectorcomponent_map.resizeArray(0);

   if (SAMRAIVectorReal<DIM,TYPE>::s_instance_counter == 0) {
      if ( !((SAMRAIVectorReal<DIM,TYPE>::s_cell_ops).isNull()) ) {
         SAMRAIVectorReal<DIM,TYPE>::s_cell_ops.setNull();
      }
      if ( !((SAMRAIVectorReal<DIM,TYPE>::s_edge_ops).isNull()) ) {
         SAMRAIVectorReal<DIM,TYPE>::s_edge_ops.setNull();
      }
      if ( !((SAMRAIVectorReal<DIM,TYPE>::s_face_ops).isNull()) ) {
         SAMRAIVectorReal<DIM,TYPE>::s_face_ops.setNull();
      }
      if ( !((SAMRAIVectorReal<DIM,TYPE>::s_node_ops).isNull()) ) {
         SAMRAIVectorReal<DIM,TYPE>::s_node_ops.setNull();
      }
      if ( !((SAMRAIVectorReal<DIM,TYPE>::s_side_ops).isNull()) ) {
         SAMRAIVectorReal<DIM,TYPE>::s_side_ops.setNull();
      }
   }

}

/*
*************************************************************************
*                                                                       *
* The following are private and cannot be used, but they are defined    *
* here for compilers that require that every template declaration have  *
* a definition (a stupid requirement, if you ask me).                   *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
SAMRAIVectorReal<DIM,TYPE>::SAMRAIVectorReal(
   const SAMRAIVectorReal<DIM,TYPE>& foo)
{
   (void) foo;  // private and not used (but included for stupid compilers)
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::operator=(
   const SAMRAIVectorReal<DIM,TYPE>& foo)
{
   (void) foo;  // private and not used (but included for stupid compilers)
}

/*
*************************************************************************
*                                                                       *
* Set name string identifier for this vector object.                    *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::setName(const std::string &name)
{
   d_vector_name = name;
}

/*
*************************************************************************
*                                                                       *
* Reset vector levels and data operations.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::resetLevels(const int coarsest_level,
                                               const int finest_level)
{
   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;
}

/*
*************************************************************************
*                                                                       *
* Create new vector with same structure as this and return new vector.  *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> >
SAMRAIVectorReal<DIM,TYPE>::cloneVector(const std::string &name) const
{

   std::string new_name = (name.empty() ? d_vector_name : name);
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > new_vec = 
      new SAMRAIVectorReal<DIM,TYPE>(new_name, 
                                       d_hierarchy, 
                                       d_coarsest_level, 
                                       d_finest_level);

   new_vec->setNumberOfComponents(d_number_components); 

   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   for (int i = 0; i < d_number_components; i++) {

      int new_id = 
         var_db->registerClonedPatchDataIndex(d_component_variable[i],
                                              d_component_data_id[i]); 

      new_vec->setComponent(i,
                            d_component_variable[i],
                            new_id,
                            d_control_volume_data_id[i], 
                            d_component_operations[i]);
   }

   return (new_vec);
}


/*
*************************************************************************
*                                                                       *
* Deallocate vector data storage and remove data indices from database. *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::freeVectorComponents()
{
   // deallocate storage for vector components
   deallocateVectorData();

   hier::VariableDatabase<DIM>* var_db = hier::VariableDatabase<DIM>::getDatabase();

   // free entries from variable database and return 
   // patch descriptor indices
   for (int i = 0; i < d_number_components; i++) {
      var_db->removePatchDataIndex(d_component_data_id[i]);
   }

   // reset variable state
   d_number_components = 0;

   d_component_variable.resizeArray(0);
   d_component_data_id.resizeArray(0);
   d_component_operations.resizeArray(0);
   d_control_volume_data_id.resizeArray(0);

   d_variableid_2_vectorcomponent_map.resizeArray(0);
}


/*
*************************************************************************
*                                                                       *
* Add new component to vector structure given a variable and the        *
* patch descriptor indices for its data and an appropriate              *
* control volume.                                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::addComponent(
   const tbox::Pointer< hier::Variable<DIM> >& var, 
   const int comp_data_id, 
   const int comp_vol_id,
   const tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> > vop) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   hier::VariableDatabase<DIM>* var_db = 
      hier::VariableDatabase<DIM>::getDatabase();
   tbox::Pointer< hier::PatchDescriptor<DIM> > patch_descriptor =
      var_db->getPatchDescriptor();
   if (!var_db->checkVariablePatchDataIndexType(var, comp_data_id)) {
      TBOX_ERROR("Error in SAMRAIVectorReal<DIM>::addComponent : "
                 << "Vector name = " << d_vector_name
                 << "\nVariable " << var->getName()
                 << " type does not match data type associated with" 
                 << " comp_data_id patch data index function argument"
                 << "\n\t var type = " << typeid(*var).name()
                 << "\n\t comp_data_id type = " 
                 << typeid(*(patch_descriptor->getPatchDataFactory(comp_data_id))).name()
                 << std::endl);
   }
    
   if (comp_vol_id >= 0) {
      if (!var_db->checkVariablePatchDataIndexType(var, comp_vol_id)) {
         TBOX_ERROR("Error in SAMRAIVectorReal<DIM>::addComponent : " 
                    << "Vector name = " << d_vector_name
                    << "\nVariable " << var->getName() 
                    << " type does not match data type associated with"
                    << " comp_vol_id patch data index function argument"
                    << "\n\t var type = " << typeid(*var).name()
                    << "\n\t comp_vol_id type = "
                    << typeid(*(patch_descriptor->getPatchDataFactory(comp_vol_id))).name()
                    << std::endl);
      }
   }
#endif

   d_number_components++;

   d_component_variable.resizeArray(d_number_components);
   d_component_data_id.resizeArray(d_number_components);
   d_component_operations.resizeArray(d_number_components);
   d_control_volume_data_id.resizeArray(d_number_components);

   hier::VariableDatabase<DIM>::getDatabase()->registerPatchDataIndex(var,
                                                                      comp_data_id);

   setComponent(d_number_components-1,
                var, 
                comp_data_id, 
                comp_vol_id,
                vop);
}

/*
*************************************************************************
*                                                                       *
* Routines to allocate and deallocate data for all vector components.   *
*                                                                       *
*************************************************************************
*/
template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::allocateVectorData(
   const double timestamp, 
   tbox::Pointer<tbox::Arena> pool)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->
                                             getPatchLevel(ln); 
      for (int i = 0; i < d_number_components; i++) {
         level->allocatePatchData(d_component_data_id[i], timestamp, pool);
      }
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::deallocateVectorData()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!d_hierarchy.isNull());
   TBOX_ASSERT(   (d_coarsest_level >= 0)
          && (d_finest_level >= d_coarsest_level)
          && (d_finest_level <= d_hierarchy->getFinestLevelNumber()) );
#endif

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      tbox::Pointer< hier::PatchLevel<DIM> > level = d_hierarchy->
                                             getPatchLevel(ln);
      for (int i = 0; i < d_number_components; i++) {
         level->deallocatePatchData(d_component_data_id[i]);
      }
   }
}


/*
*************************************************************************
*                                                                       *
* Print Vector attributes and data.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::print(std::ostream& s, bool interior_only) const
{
   s << "\nVector : " << getName() << std::endl;
   s << "coarsest level = " << d_coarsest_level 
     << " : finest level = " << d_finest_level << std::endl;
   s << "d_number_components = " << d_number_components << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Printing data components on level " << ln << std::endl;
      for (int i = 0; i < d_number_components; i++) {
         s << "Vector component index = " << i << std::endl; 
	 d_component_operations[i]->resetLevels(ln,ln);
         d_component_operations[i]->printData(d_component_data_id[i],
					      s,
					      interior_only);
      }
   }
   return;
}

/*
*************************************************************************
*                                                                       *
* Private member functions to set the number of vector components       *
* to set individual components.   These routines are used when cloning  *
* vectors and/or adding components.                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::setNumberOfComponents(int num_comp)
{
   d_number_components = num_comp;

   d_component_variable.resizeArray(d_number_components);
   d_component_data_id.resizeArray(d_number_components);
   d_component_operations.resizeArray(d_number_components);
   d_control_volume_data_id.resizeArray(d_number_components);
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::setComponent(
   const int comp_id,
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const int data_id,
   const int vol_id,
   const tbox::Pointer< math::HierarchyDataOpsReal<DIM,TYPE> > vop)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id < d_number_components);
#endif

   d_component_variable[comp_id] = var;
   d_component_data_id[comp_id]  = data_id;
   if ( vop.isNull() ) {

      const tbox::Pointer< pdat::CellVariable<DIM,TYPE> > cellvar(var);
      const tbox::Pointer< pdat::EdgeVariable<DIM,TYPE> > edgevar(var);
      const tbox::Pointer< pdat::FaceVariable<DIM,TYPE> > facevar(var);
      const tbox::Pointer< pdat::NodeVariable<DIM,TYPE> > nodevar(var);
      const tbox::Pointer< pdat::SideVariable<DIM,TYPE> > sidevar(var);

      if ( !(cellvar.isNull()) ) {
         if (!SAMRAIVectorReal<DIM,TYPE>::s_cell_ops) {
            SAMRAIVectorReal<DIM,TYPE>::s_cell_ops =
            new math::HierarchyCellDataOpsReal<DIM,TYPE>(d_hierarchy,
                                                     d_coarsest_level,
                                                     d_finest_level);
         }
         d_component_operations[comp_id] = 
            SAMRAIVectorReal<DIM,TYPE>::s_cell_ops;
      } else if ( !(edgevar.isNull()) ) {
         if (!SAMRAIVectorReal<DIM,TYPE>::s_edge_ops) {
            SAMRAIVectorReal<DIM,TYPE>::s_edge_ops =
            new math::HierarchyEdgeDataOpsReal<DIM,TYPE>(d_hierarchy,
                                                     d_coarsest_level,
                                                     d_finest_level);
         }
         d_component_operations[comp_id] = 
            SAMRAIVectorReal<DIM,TYPE>::s_edge_ops;
      } else if ( !(facevar.isNull()) ) {
         if (!SAMRAIVectorReal<DIM,TYPE>::s_face_ops) {
            SAMRAIVectorReal<DIM,TYPE>::s_face_ops =
            new math::HierarchyFaceDataOpsReal<DIM,TYPE>(d_hierarchy,
                                                     d_coarsest_level,
                                                     d_finest_level);
         }
         d_component_operations[comp_id] = 
            SAMRAIVectorReal<DIM,TYPE>::s_face_ops;
      } else if ( !(nodevar.isNull()) ) {
         if (!SAMRAIVectorReal<DIM,TYPE>::s_node_ops) {
            SAMRAIVectorReal<DIM,TYPE>::s_node_ops =
            new math::HierarchyNodeDataOpsReal<DIM,TYPE>(d_hierarchy,
                                                     d_coarsest_level,
                                                     d_finest_level);
         }
         d_component_operations[comp_id] = 
            SAMRAIVectorReal<DIM,TYPE>::s_node_ops;
      } else if ( !(sidevar.isNull()) ) {
         if (!SAMRAIVectorReal<DIM,TYPE>::s_side_ops) {
            SAMRAIVectorReal<DIM,TYPE>::s_side_ops =
            new math::HierarchySideDataOpsReal<DIM,TYPE>(d_hierarchy,
                                                     d_coarsest_level,
                                                     d_finest_level);
         }
         d_component_operations[comp_id] = 
            SAMRAIVectorReal<DIM,TYPE>::s_side_ops;
      }
   } else {
      d_component_operations[comp_id]  = vop;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!(d_component_operations[comp_id].isNull()));
#endif

   d_control_volume_data_id[comp_id] = vol_id;

   int var_id = var->getInstanceIdentifier();

   int oldsize = d_variableid_2_vectorcomponent_map.getSize();
   int newsize = var_id + 1;
   if (oldsize < newsize) {
      newsize = tbox::MathUtilities<int>::Max(
                oldsize + DESCRIPTOR_ID_ARRAY_SCRATCH_SPACE, newsize );
      d_variableid_2_vectorcomponent_map.resizeArray(newsize);
      for (int i = oldsize; i < newsize; i++) {
         d_variableid_2_vectorcomponent_map[i] = -1;
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_variableid_2_vectorcomponent_map[var_id] == -1);
#endif

   d_variableid_2_vectorcomponent_map[var_id] = comp_id;
}


/*
*************************************************************************
*                                                                       *
* The remaining functions are basic vector kernel routines.             *
* The operation for each component is performed by its hierarchy        *
* data operation object.                                                *
*                                                                       *
*************************************************************************
*/

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::copyVector(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > src_vec,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->
      copyData(d_component_data_id[i],
               src_vec->getComponentDescriptorIndex(i),
	       interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::swapVectors(
   tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > other)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->
      swapData(d_component_data_id[i],
               other->getComponentDescriptorIndex(i));
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::setToScalar(const TYPE& alpha,
					       const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->setToScalar(d_component_data_id[i],
                                             alpha,
					     interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::scale(
   const TYPE& alpha, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->scale(d_component_data_id[i],
                                       alpha,
                                       x->getComponentDescriptorIndex(i),
				       interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::addScalar(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x, 
   const TYPE& alpha,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->addScalar(d_component_data_id[i],
                                           x->getComponentDescriptorIndex(i),
                                           alpha,
					   interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::add(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->add(d_component_data_id[i],
                                     x->getComponentDescriptorIndex(i),
                                     y->getComponentDescriptorIndex(i),
				     interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::subtract(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->subtract(d_component_data_id[i],
                                          x->getComponentDescriptorIndex(i),
                                          y->getComponentDescriptorIndex(i),
					  interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::multiply(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->multiply(d_component_data_id[i],
                                          x->getComponentDescriptorIndex(i),
                                          y->getComponentDescriptorIndex(i),
					  interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::divide(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->divide(d_component_data_id[i],
                                        x->getComponentDescriptorIndex(i),
                                        y->getComponentDescriptorIndex(i),
					interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::reciprocal(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->reciprocal(d_component_data_id[i],
                                            x->getComponentDescriptorIndex(i),
					    interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::linearSum(
   const TYPE& alpha, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   const TYPE& beta, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > y,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->linearSum(d_component_data_id[i],
                                           alpha,
                                           x->getComponentDescriptorIndex(i),
                                           beta,
                                           y->getComponentDescriptorIndex(i),
					   interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::axpy(
   const TYPE& alpha, 
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > y,
   const bool interior_only) 
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->axpy(d_component_data_id[i],
                                      alpha,
                                      x->getComponentDescriptorIndex(i),
                                      y->getComponentDescriptorIndex(i),
				      interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::abs(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->abs(d_component_data_id[i],
                                     x->getComponentDescriptorIndex(i),
				     interior_only);
   }
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::setRandomValues(const TYPE& width,
                                                   const TYPE& low,
						   const bool interior_only)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->setRandomValues(d_component_data_id[i],
                                                 width,
                                                 low,
						 interior_only);
   }
}

template<int DIM, class TYPE>
double SAMRAIVectorReal<DIM,TYPE>::L1Norm(bool local_only) const
{
   double norm = 0.0;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      norm += d_component_operations[i]->L1Norm(d_component_data_id[i],
                                                d_control_volume_data_id[i]);
   }

   return( norm );
}

template<int DIM, class TYPE>
double SAMRAIVectorReal<DIM,TYPE>::L2Norm(bool local_only) const
{
   double norm_squared = 0.0;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      double comp_norm = 
         d_component_operations[i]->L2Norm(d_component_data_id[i],
                                           d_control_volume_data_id[i]);
      norm_squared += comp_norm * comp_norm;
   }

   return( sqrt(norm_squared) );
}

template<int DIM, class TYPE>
double SAMRAIVectorReal<DIM,TYPE>::weightedL2Norm(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > wgt) const
{
   double norm_squared = 0.0;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      double comp_norm = d_component_operations[i]->weightedL2Norm(
                         d_component_data_id[i],
                         wgt->getComponentDescriptorIndex(i),
                         d_control_volume_data_id[i]);
      norm_squared += comp_norm * comp_norm;
   }

   return( sqrt(norm_squared) );
}

template<int DIM, class TYPE>
double SAMRAIVectorReal<DIM,TYPE>::RMSNorm() const
{
   double num = L2Norm();

   double denom = 0.0;
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      if ( d_control_volume_data_id[i] < 0 ) {
         denom += double( d_component_operations[i]->
                          numberOfEntries(d_component_data_id[i], true) );
      } else {
         denom += d_component_operations[i]->
                  sumControlVolumes(d_component_data_id[i],
                                    d_control_volume_data_id[i]);
      }
   }
  
   double norm = 0.0; 
   if (denom > 0.0) norm = num/sqrt(denom);
   return( norm );
}

template<int DIM, class TYPE>
double SAMRAIVectorReal<DIM,TYPE>::weightedRMSNorm(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > wgt) const
{
   double num = weightedL2Norm(wgt);

   double denom = 0.0;
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      if ( d_control_volume_data_id[i] < 0 ) {
         denom += double( d_component_operations[i]->
                          numberOfEntries(d_component_data_id[i], true) );
      } else {
         denom += d_component_operations[i]->
                  sumControlVolumes(d_component_data_id[i],
                                    d_control_volume_data_id[i]);
      }
   }

   double norm = 0.0;
   if (denom > 0.0) norm = num/sqrt(denom); 
   return( norm );
}

template<int DIM, class TYPE>
double SAMRAIVectorReal<DIM,TYPE>::maxNorm(bool local_only) const
{
   double norm = 0.0;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      norm = tbox::MathUtilities<double>::Max( norm, 
                   d_component_operations[i]->maxNorm(
                                              d_component_data_id[i], 
                                              d_control_volume_data_id[i]) );
   }

   return( norm );
}

template<int DIM, class TYPE>
TYPE SAMRAIVectorReal<DIM,TYPE>::dot(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   bool local_only) const 
{
   TYPE dprod = 0.0;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      dprod += d_component_operations[i]->dot(d_component_data_id[i],
                                              x->getComponentDescriptorIndex(i),
                                              d_control_volume_data_id[i]);
   }

   return( dprod );
}

template<int DIM, class TYPE>
int SAMRAIVectorReal<DIM,TYPE>::computeConstrProdPos(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x) const
{
   int test = 1;

   int i = 0;
   while (test && (i < d_number_components)) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      test = tbox::MathUtilities<int>::Min( test,
                   d_component_operations[i]->
                   computeConstrProdPos(d_component_data_id[i],
                                        x->getComponentDescriptorIndex(i),
                                        d_control_volume_data_id[i]) );
      i++;
   }

   return( test );
}

template<int DIM, class TYPE>
void SAMRAIVectorReal<DIM,TYPE>::compareToScalar(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x,
   const TYPE& alpha)
{
   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      d_component_operations[i]->
      compareToScalar(d_component_data_id[i],
                      x->getComponentDescriptorIndex(i),
                      alpha,
                      d_control_volume_data_id[i]);
   }
}

template<int DIM, class TYPE>
int SAMRAIVectorReal<DIM,TYPE>::testReciprocal(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > x)
{
   int test = 1;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      test = tbox::MathUtilities<int>::Min( test,
                   d_component_operations[i]->
                   testReciprocal(d_component_data_id[i],
                                  x->getComponentDescriptorIndex(i),
                                  d_control_volume_data_id[i]) );
   }

   return( test );
}

template<int DIM, class TYPE>
TYPE SAMRAIVectorReal<DIM,TYPE>::maxPointwiseDivide(
   const tbox::Pointer< SAMRAIVectorReal<DIM,TYPE> > denom) const
{
   TYPE max = 0.0;

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      TYPE component_max =
	 d_component_operations[i]->maxPointwiseDivide(d_component_data_id[i],
         denom->getComponentDescriptorIndex(i),
						       true);
      max = tbox::MathUtilities<double>::Max( max, component_max );
   }

   max = tbox::SAMRAI_MPI::maxReduction(max);
   return( max );
}


template<int DIM, class TYPE>
TYPE SAMRAIVectorReal<DIM,TYPE>::min(const bool interior_only) const
{
   TYPE minval = tbox::MathUtilities<TYPE>::getMax();

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      minval = tbox::MathUtilities<TYPE>::Min(
               minval,
               d_component_operations[i]->min(d_component_data_id[i],
					      interior_only) );
   }

   return( minval );
}

template<int DIM, class TYPE>
TYPE SAMRAIVectorReal<DIM,TYPE>::max(const bool interior_only) const
{
   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();

   for (int i = 0; i < d_number_components; i++) {
      d_component_operations[i]->resetLevels(d_coarsest_level, d_finest_level);
      maxval = tbox::MathUtilities<TYPE>::Max(
               maxval,
               d_component_operations[i]->max(d_component_data_id[i],
					      interior_only) );
   }

   return( maxval );
}


}
}
#endif
