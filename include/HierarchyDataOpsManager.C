//
// File:        HierarchyDataOpsManager.C
// Package:     SAMRAI mathops
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Singleton manager for hierarchy data operation objects.
//

#ifndef included_math_HierarchyDataOpsManager_C
#define included_math_HierarchyDataOpsManager_C

#include "HierarchyDataOpsManager.h"

#include "HierarchyCellDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"

#ifdef HAVE_DCOMPLEX
#include "HierarchyCellDataOpsComplex.h"
#include "HierarchyFaceDataOpsComplex.h"
#include "HierarchyNodeDataOpsComplex.h"
#include "HierarchySideDataOpsComplex.h"
#include "HierarchyEdgeDataOpsComplex.h"
#endif

#include "HierarchyCellDataOpsInteger.h"
#include "HierarchyFaceDataOpsInteger.h"
#include "HierarchyNodeDataOpsInteger.h"
#include "HierarchySideDataOpsInteger.h"
#include "HierarchyEdgeDataOpsInteger.h"

#include "CellVariable.h"
#include "FaceVariable.h"
#include "NodeVariable.h"
#include "SideVariable.h"
#include "EdgeVariable.h"
#include "tbox/Complex.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace math {

/*
*************************************************************************
*                                                                       *
* Static members for Singleton hierarchy operation manager class.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> HierarchyDataOpsManager<DIM>* 
HierarchyDataOpsManager<DIM>::s_pdat_op_manager_instance = 
                           ((HierarchyDataOpsManager<DIM>*) NULL);
template<int DIM> bool HierarchyDataOpsManager<DIM>::s_registered_callback = false;

template<int DIM> HierarchyDataOpsManager<DIM>* HierarchyDataOpsManager<DIM>::getManager()
{
   if (!s_pdat_op_manager_instance) {
      s_pdat_op_manager_instance = new HierarchyDataOpsManager<DIM>();
   }
   if (!s_registered_callback) {
      tbox::ShutdownRegistry::registerShutdownRoutine(freeManager,
                       tbox::ShutdownRegistry::priorityHierarchyDataOpsManager);
      s_registered_callback = true;
   }
   return(s_pdat_op_manager_instance);
}

template<int DIM> void HierarchyDataOpsManager<DIM>::freeManager()
{
   if (s_pdat_op_manager_instance) delete s_pdat_op_manager_instance;
   s_pdat_op_manager_instance = ((HierarchyDataOpsManager<DIM>*) NULL);
}

template<int DIM> void HierarchyDataOpsManager<DIM>::registerSingletonSubclassInstance(
   HierarchyDataOpsManager<DIM>* subclass_instance)
{
   if (!s_pdat_op_manager_instance) {
      s_pdat_op_manager_instance = subclass_instance;
      if (!s_registered_callback) {
         tbox::ShutdownRegistry::registerShutdownRoutine(freeManager,
                       tbox::ShutdownRegistry::priorityHierarchyDataOpsManager);
         s_registered_callback = true;
      }
   } else {
      TBOX_ERROR("HierarchyDataOpsManager<DIM> internal error...\n"
                 << "Attempting to set Singleton instance to subclass instance,"
                 << "\n but Singleton instance already set." << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Empty constructor and destructor for hierarchy operation manager.     * 
*                                                                       *
*************************************************************************
*/

template<int DIM>  HierarchyDataOpsManager<DIM>::HierarchyDataOpsManager() 
{
}
   
template<int DIM>  HierarchyDataOpsManager<DIM>::~HierarchyDataOpsManager()
{
}

/*!
  Return pointer to operation object for a double variable
  on the given hierarchy.

  If a unique operator object is not requested, and if one already
  exists for the hierarchy and variable specified, the existing one
  will be created and returned.  Otherwise, a new one is created.
  Objects created created for unique requests will not be used later
  when an equivalent request is made.
*/

template<int DIM> tbox::Pointer< HierarchyDataOpsReal<DIM,double> >
HierarchyDataOpsManager<DIM>::getOperationsDouble(
   const tbox::Pointer< hier::Variable<DIM> >& variable,
   tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
   bool get_unique )
{
  const tbox::Pointer< pdat::CellVariable<DIM,double> > cellvar(variable);
  const tbox::Pointer< pdat::FaceVariable<DIM,double> > facevar(variable);
  const tbox::Pointer< pdat::NodeVariable<DIM,double> > nodevar(variable);
  const tbox::Pointer< pdat::SideVariable<DIM,double> > sidevar(variable);
  const tbox::Pointer< pdat::EdgeVariable<DIM,double> > edgevar(variable);

  tbox::Pointer< HierarchyDataOpsReal<DIM,double> > ops;

  if ( !(cellvar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyCellDataOpsReal<DIM,double>(hierarchy);
    }
    else {
      const int n = d_cell_ops_double.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_cell_ops_double[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_cell_ops_double[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyCellDataOpsReal<DIM,double>(hierarchy);
	d_cell_ops_double.resizeArray(n+1);
	d_cell_ops_double[n] = ops;
      }
    }

  } else if ( !(facevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyFaceDataOpsReal<DIM,double>(hierarchy);
    }
    else {
      const int n = d_face_ops_double.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_face_ops_double[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_face_ops_double[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyFaceDataOpsReal<DIM,double>(hierarchy);
	d_face_ops_double.resizeArray(n+1);
	d_face_ops_double[n] = ops;
      }
    }

  } else if ( !(nodevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyNodeDataOpsReal<DIM,double>(hierarchy);
    }
    else {
      const int n = d_node_ops_double.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_node_ops_double[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_node_ops_double[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyNodeDataOpsReal<DIM,double>(hierarchy);
	d_node_ops_double.resizeArray(n+1);
	d_node_ops_double[n] = ops;
      }
    }

  } else if ( !(sidevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchySideDataOpsReal<DIM,double>(hierarchy);
    }
    else {
      const int n = d_side_ops_double.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_side_ops_double[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_side_ops_double[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchySideDataOpsReal<DIM,double>(hierarchy);
	d_side_ops_double.resizeArray(n+1);
	d_side_ops_double[n] = ops;
      }
    }

  } else if ( !(edgevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyEdgeDataOpsReal<DIM,double>(hierarchy);
    }
    else {
      const int n = d_edge_ops_double.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_edge_ops_double[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_edge_ops_double[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyEdgeDataOpsReal<DIM,double>(hierarchy);
	d_edge_ops_double.resizeArray(n+1);
	d_edge_ops_double[n] = ops;
      }
    }

  }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager<DIM> internal error...\n"
                 << "Operations for variable " << variable->getName() 
                 << " not defined.");
   }

   return( ops ); 
}

/*!
  Return pointer to operation object for a float variable
  on the given hierarchy.

  If a unique operator object is not requested, and if one already
  exists for the hierarchy and variable specified, the existing one
  will be created and returned.  Otherwise, a new one is created.
  Objects created created for unique requests will not be used later
  when an equivalent request is made.
*/

#ifdef HAVE_FLOAT
template<int DIM> tbox::Pointer< HierarchyDataOpsReal<DIM,float> >
HierarchyDataOpsManager<DIM>::getOperationsFloat(
   const tbox::Pointer< hier::Variable<DIM> >& variable,
   tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
   bool get_unique )
{
  const tbox::Pointer< pdat::CellVariable<DIM,float> > cellvar(variable);
  const tbox::Pointer< pdat::FaceVariable<DIM,float> > facevar(variable);
  const tbox::Pointer< pdat::NodeVariable<DIM,float> > nodevar(variable);
  const tbox::Pointer< pdat::SideVariable<DIM,float> > sidevar(variable);
  const tbox::Pointer< pdat::EdgeVariable<DIM,float> > edgevar(variable);

  tbox::Pointer< HierarchyDataOpsReal<DIM,float> > ops;

  if ( !(cellvar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyCellDataOpsReal<DIM,float>(hierarchy);
    }
    else {
      const int n = d_cell_ops_float.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_cell_ops_float[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_cell_ops_float[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyCellDataOpsReal<DIM,float>(hierarchy);
	d_cell_ops_float.resizeArray(n+1);
	d_cell_ops_float[n] = ops;
      }
    }

  } else if ( !(facevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyFaceDataOpsReal<DIM,float>(hierarchy);
    }
    else {
      const int n = d_face_ops_float.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_face_ops_float[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_face_ops_float[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyFaceDataOpsReal<DIM,float>(hierarchy);
	d_face_ops_float.resizeArray(n+1);
	d_face_ops_float[n] = ops;
      }
    }

  } else if ( !(nodevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyNodeDataOpsReal<DIM,float>(hierarchy);
    }
    else {
      const int n = d_node_ops_float.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_node_ops_float[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_node_ops_float[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyNodeDataOpsReal<DIM,float>(hierarchy);
	d_node_ops_float.resizeArray(n+1);
	d_node_ops_float[n] = ops;
      }
    }

  } else if ( !(sidevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchySideDataOpsReal<DIM,float>(hierarchy);
    }
    else {
      const int n = d_side_ops_float.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_side_ops_float[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_side_ops_float[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchySideDataOpsReal<DIM,float>(hierarchy);
	d_side_ops_float.resizeArray(n+1);
	d_side_ops_float[n] = ops;
      }
    }

  } else if ( !(edgevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyEdgeDataOpsReal<DIM,float>(hierarchy);
    }
    else {
      const int n = d_edge_ops_float.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_edge_ops_float[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_edge_ops_float[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyEdgeDataOpsReal<DIM,float>(hierarchy);
	d_edge_ops_float.resizeArray(n+1);
	d_edge_ops_float[n] = ops;
      }
    }

  }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager<DIM> internal error...\n"
                 << "Operations for variable " << variable->getName() 
                 << " not defined.");
   }

   return( ops ); 
}
#endif

/*!
  Return pointer to operation object for a complex variable
  on the given hierarchy.

  If a unique operator object is not requested, and if one already
  exists for the hierarchy and variable specified, the existing one
  will be created and returned.  Otherwise, a new one is created.
  Objects created created for unique requests will not be used later
  when an equivalent request is made.
*/

#ifdef HAVE_DCOMPLEX
template<int DIM> tbox::Pointer< HierarchyDataOpsComplex<DIM> >
HierarchyDataOpsManager<DIM>::getOperationsComplex(
   const tbox::Pointer< hier::Variable<DIM> >& variable,
   tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
   bool get_unique )
{
  const tbox::Pointer< pdat::CellVariable<DIM,dcomplex> > cellvar(variable);
  const tbox::Pointer< pdat::FaceVariable<DIM,dcomplex> > facevar(variable);
  const tbox::Pointer< pdat::NodeVariable<DIM,dcomplex> > nodevar(variable);
  const tbox::Pointer< pdat::SideVariable<DIM,dcomplex> > sidevar(variable);
  const tbox::Pointer< pdat::EdgeVariable<DIM,dcomplex> > edgevar(variable);

  tbox::Pointer< HierarchyDataOpsComplex<DIM> > ops;

  if ( !(cellvar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyCellDataOpsComplex<DIM>(hierarchy);
    }
    else {
      const int n = d_cell_ops_double.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_cell_ops_complex[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_cell_ops_complex[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyCellDataOpsComplex<DIM>(hierarchy);
	d_cell_ops_complex.resizeArray(n+1);
	d_cell_ops_complex[n] = ops;
      }
    }

  } else if ( !(facevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyFaceDataOpsComplex<DIM>(hierarchy);
    }
    else {
      const int n = d_face_ops_complex.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_face_ops_complex[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_face_ops_complex[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyFaceDataOpsComplex<DIM>(hierarchy);
	d_face_ops_complex.resizeArray(n+1);
	d_face_ops_complex[n] = ops;
      }
    }

  } else if ( !(nodevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyNodeDataOpsComplex<DIM>(hierarchy);
    }
    else {
      const int n = d_node_ops_complex.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_node_ops_complex[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_node_ops_complex[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyNodeDataOpsComplex<DIM>(hierarchy);
	d_node_ops_complex.resizeArray(n+1);
	d_node_ops_complex[n] = ops;
      }
    }

  } else if ( !(sidevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchySideDataOpsComplex<DIM>(hierarchy);
    }
    else {
      const int n = d_side_ops_complex.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_side_ops_complex[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_side_ops_complex[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchySideDataOpsComplex<DIM>(hierarchy);
	d_side_ops_complex.resizeArray(n+1);
	d_side_ops_complex[n] = ops;
      }
    }

  } else if ( !(edgevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyEdgeDataOpsComplex<DIM>(hierarchy);
    }
    else {
      const int n = d_edge_ops_complex.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_edge_ops_complex[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_edge_ops_complex[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyEdgeDataOpsComplex<DIM>(hierarchy);
	d_edge_ops_complex.resizeArray(n+1);
	d_edge_ops_complex[n] = ops;
      }
    }

  }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager<DIM> internal error...\n"
                 << "Operations for variable " << variable->getName() 
                 << " not defined.");
   }

   return( ops ); 
}
#endif

/*!
  Return pointer to operation object for an integer variable
  on the given hierarchy.

  If a unique operator object is not requested, and if one already
  exists for the hierarchy and variable specified, the existing one
  will be created and returned.  Otherwise, a new one is created.
  Objects created created for unique requests will not be used later
  when an equivalent request is made.
*/

template<int DIM> tbox::Pointer< HierarchyDataOpsInteger<DIM> >
HierarchyDataOpsManager<DIM>::getOperationsInteger(
   const tbox::Pointer< hier::Variable<DIM> >& variable,
   tbox::Pointer< hier::PatchHierarchy<DIM> >& hierarchy,
   bool get_unique )
{
  const tbox::Pointer< pdat::CellVariable<DIM,int> > cellvar(variable);
  const tbox::Pointer< pdat::FaceVariable<DIM,int> > facevar(variable);
  const tbox::Pointer< pdat::NodeVariable<DIM,int> > nodevar(variable);
  const tbox::Pointer< pdat::SideVariable<DIM,int> > sidevar(variable);
  const tbox::Pointer< pdat::EdgeVariable<DIM,int> > edgevar(variable);

  tbox::Pointer< HierarchyDataOpsInteger<DIM> > ops;

  if ( !(cellvar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyCellDataOpsInteger<DIM>(hierarchy);
    }
    else {
      const int n = d_cell_ops_int.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_cell_ops_int[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_cell_ops_int[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyCellDataOpsInteger<DIM>(hierarchy);
	d_cell_ops_int.resizeArray(n+1);
	d_cell_ops_int[n] = ops;
      }
    }

  } else if ( !(facevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyFaceDataOpsInteger<DIM>(hierarchy);
    }
    else {
      const int n = d_face_ops_int.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_face_ops_int[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_face_ops_int[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyFaceDataOpsInteger<DIM>(hierarchy);
	d_face_ops_int.resizeArray(n+1);
	d_face_ops_int[n] = ops;
      }
    }

  } else if ( !(nodevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyNodeDataOpsInteger<DIM>(hierarchy);
    }
    else {
      const int n = d_node_ops_int.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_node_ops_int[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_node_ops_int[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyNodeDataOpsInteger<DIM>(hierarchy);
	d_node_ops_int.resizeArray(n+1);
	d_node_ops_int[n] = ops;
      }
    }

  } else if ( !(sidevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchySideDataOpsInteger<DIM>(hierarchy);
    }
    else {
      const int n = d_side_ops_int.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_side_ops_int[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_side_ops_int[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchySideDataOpsInteger<DIM>(hierarchy);
	d_side_ops_int.resizeArray(n+1);
	d_side_ops_int[n] = ops;
      }
    }

  } else if ( !(edgevar.isNull()) ) {

    if ( get_unique ) {
      ops = new HierarchyEdgeDataOpsInteger<DIM>(hierarchy);
    }
    else {
      const int n = d_edge_ops_int.getSize();
      for ( int i=0; i<n && ops.isNull(); ++i ) {
        if ( hierarchy != d_edge_ops_int[i]->getPatchHierarchy() ) continue;
	// A compatible operator has been found at i.
	ops = d_edge_ops_int[i];
      }
      if (!ops) {
	// No compatible operator has been found.
	ops = new HierarchyEdgeDataOpsInteger<DIM>(hierarchy);
	d_edge_ops_int.resizeArray(n+1);
	d_edge_ops_int[n] = ops;
      }
    }

  }

   if (!ops) {
      TBOX_ERROR("HierarchyDataOpsManager<DIM> internal error...\n"
                 << "Operations for variable " << variable->getName() 
                 << " not defined.");
   }

   return( ops ); 
}


}
}
#endif
