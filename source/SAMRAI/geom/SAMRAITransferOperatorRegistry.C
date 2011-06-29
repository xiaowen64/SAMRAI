/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Singleton registry for all tranfer operators. 
 *
 ************************************************************************/

#ifndef included_hier_SAMRAITransferOperatorRegistry_C
#define included_hier_SAMRAITransferOperatorRegistry_C

#include "SAMRAI/geom/SAMRAITransferOperatorRegistry.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/pdat/NodeComplexInjection.h"
#include "SAMRAI/pdat/NodeDoubleInjection.h"
#include "SAMRAI/pdat/NodeFloatInjection.h"
#include "SAMRAI/pdat/NodeIntegerInjection.h"
#include "SAMRAI/pdat/OuternodeDoubleConstantCoarsen.h"
#include "SAMRAI/pdat/CellComplexConstantRefine.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellFloatConstantRefine.h"
#include "SAMRAI/pdat/CellIntegerConstantRefine.h"
#include "SAMRAI/pdat/EdgeComplexConstantRefine.h"
#include "SAMRAI/pdat/EdgeDoubleConstantRefine.h"
#include "SAMRAI/pdat/EdgeFloatConstantRefine.h"
#include "SAMRAI/pdat/EdgeIntegerConstantRefine.h"
#include "SAMRAI/pdat/FaceComplexConstantRefine.h"
#include "SAMRAI/pdat/FaceDoubleConstantRefine.h"
#include "SAMRAI/pdat/FaceFloatConstantRefine.h"
#include "SAMRAI/pdat/FaceIntegerConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceComplexConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceDoubleConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceFloatConstantRefine.h"
#include "SAMRAI/pdat/OuterfaceIntegerConstantRefine.h"
#include "SAMRAI/pdat/SideComplexConstantRefine.h"
#include "SAMRAI/pdat/SideDoubleConstantRefine.h"
#include "SAMRAI/pdat/SideFloatConstantRefine.h"
#include "SAMRAI/pdat/SideIntegerConstantRefine.h"
#include "SAMRAI/pdat/CellComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/CellDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/CellFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/EdgeFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/FaceFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/NodeFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OuterfaceFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/OutersideFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideComplexLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideDoubleLinearTimeInterpolateOp.h"
#include "SAMRAI/pdat/SideFloatLinearTimeInterpolateOp.h"
#include "SAMRAI/geom/CartesianCellComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellFloatWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeFloatWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceFloatWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceFloatWeightedAverage.h"
#include "SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideFloatWeightedAverage.h"
#include "SAMRAI/geom/SkeletonCoarsen.h"
#include "SAMRAI/geom/CartesianCellComplexConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianEdgeFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianFaceFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianSideDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianSideFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeFloatLinearRefine.h"
#include "SAMRAI/geom/SkeletonRefine.h"

namespace SAMRAI {
namespace geom {

/*
 *************************************************************************
 *                                                                       *
 * Constructor and destructor for TransferOperatorRegistry objects.      *
 *                                                                       *
 *************************************************************************
 */

SAMRAITransferOperatorRegistry::SAMRAITransferOperatorRegistry(
   const tbox::Dimension& dim):
   TransferOperatorRegistry(dim)
{
}

SAMRAITransferOperatorRegistry::~SAMRAITransferOperatorRegistry()
{
}

tbox::Pointer<hier::CoarsenOperator>
SAMRAITransferOperatorRegistry::buildCoarsenOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name)
{
   tbox::Pointer<hier::CoarsenOperator> coarsen_op(NULL);
   if (op_name == "CONSERVATIVE_COARSEN") {
      if (!tbox::Pointer<pdat::CellVariable<dcomplex> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianCellComplexWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianCellDoubleWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<float> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianCellFloatWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<dcomplex> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianEdgeComplexWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianEdgeDoubleWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<float> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianEdgeFloatWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<dcomplex> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianFaceComplexWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianFaceDoubleWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<float> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianFaceFloatWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<dcomplex> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianOuterfaceComplexWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianOuterfaceDoubleWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<float> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianOuterfaceFloatWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::OutersideVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianOutersideDoubleWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<dcomplex> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianSideComplexWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianSideDoubleWeightedAverage(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<float> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new CartesianSideFloatWeightedAverage(getDim()));
      }
      else {
         TBOX_ERROR("Unsupported variable centering for coarsen operator");
      }
   }
   else if (op_name == "SKELETON_COARSEN") {
      coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
         new SkeletonCoarsen(getDim()));
   }
   else if (op_name == "CONSTANT_COARSEN") {
      if (!tbox::Pointer<pdat::NodeVariable<dcomplex> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new pdat::NodeComplexInjection(getDim()));
      }
      else if (!tbox::Pointer<pdat::NodeVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new pdat::NodeDoubleInjection(getDim()));
      }
      else if (!tbox::Pointer<pdat::NodeVariable<float> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new pdat::NodeFloatInjection(getDim()));
      }
      else if (!tbox::Pointer<pdat::NodeVariable<int> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new pdat::NodeIntegerInjection(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuternodeVariable<double> >(var).isNull()) {
         coarsen_op = tbox::Pointer<hier::CoarsenOperator>(
            new pdat::OuternodeDoubleConstantCoarsen(getDim()));
      }
      else {
         TBOX_ERROR("Unsupported variable centering for coarsen operator");
      }
   }
   else {
      TBOX_ERROR("Unknown coarsen operator name");
   }
   addCoarsenOperator(coarsen_op);
   return coarsen_op;
}

tbox::Pointer<hier::RefineOperator>
SAMRAITransferOperatorRegistry::buildRefineOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name)
{
   tbox::Pointer<hier::RefineOperator> refine_op(NULL);
   if (op_name == "CONSERVATIVE_LINEAR_REFINE") {
      if (!tbox::Pointer<pdat::CellVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianCellComplexConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianCellDoubleConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianCellFloatConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianEdgeDoubleConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianEdgeFloatConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianFaceDoubleConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianFaceFloatConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianSideDoubleConservativeLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianSideFloatConservativeLinearRefine(getDim()));
      }
      else {
         TBOX_ERROR("Unsupported variable centering for refine operator");
      }
   }
   else if (op_name == "LINEAR_REFINE") {
      if (!tbox::Pointer<pdat::CellVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianCellComplexLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianCellDoubleLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianCellFloatLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::NodeVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianNodeComplexLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::NodeVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianNodeDoubleLinearRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::NodeVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new CartesianNodeFloatLinearRefine(getDim()));
      }
      else {
         TBOX_ERROR("Unsupported variable centering for refine operator");
      }
   }
   else if (op_name == "SKELETON_REFINE") {
      refine_op = tbox::Pointer<hier::RefineOperator>(
         new SkeletonRefine(getDim()));
   }
   else if (op_name == "CONSTANT_REFINE") {
      if (!tbox::Pointer<pdat::CellVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::CellComplexConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::CellDoubleConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::CellFloatConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::CellVariable<int> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::CellIntegerConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::EdgeComplexConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::EdgeDoubleConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::EdgeFloatConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<int> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::EdgeIntegerConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::FaceComplexConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::FaceDoubleConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::FaceFloatConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::FaceVariable<int> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::FaceIntegerConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::OuterfaceComplexConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::OuterfaceDoubleConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::OuterfaceFloatConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<int> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::OuterfaceIntegerConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<dcomplex> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::SideComplexConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<double> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::SideDoubleConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<float> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::SideFloatConstantRefine(getDim()));
      }
      else if (!tbox::Pointer<pdat::SideVariable<int> >(var).isNull()) {
         refine_op = tbox::Pointer<hier::RefineOperator>(
            new pdat::SideIntegerConstantRefine(getDim()));
      }
      else {
         TBOX_ERROR("Unsupported variable centering for refine operator");
      }
   }
   else {
      TBOX_ERROR("Unknown coarsen operator name");
   }
   addRefineOperator(refine_op);
   return refine_op;
}

tbox::Pointer<hier::TimeInterpolateOperator>
SAMRAITransferOperatorRegistry::buildTimeInterpolateOperator(
   const tbox::Pointer<hier::Variable>& var,
   const std::string& op_name)
{
   tbox::Pointer<hier::TimeInterpolateOperator> time_op(NULL);
   if (op_name == "STD_LINEAR_TIME_INTERPOLATE") {
      if (!tbox::Pointer<pdat::CellVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::CellComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::CellVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::CellDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::CellVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::CellFloatLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::EdgeComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::EdgeDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::EdgeVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::EdgeFloatLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::FaceVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::FaceComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::FaceVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::FaceDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::FaceVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::FaceFloatLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::NodeVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::NodeComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::NodeVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::NodeDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::NodeVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::NodeFloatLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::OuterfaceComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::OuterfaceDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::OuterfaceVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::OuterfaceFloatLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::OutersideVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::OutersideComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::OutersideVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::OutersideDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::OutersideVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::OutersideFloatLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::SideVariable<dcomplex> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::SideComplexLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::SideVariable<double> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::SideDoubleLinearTimeInterpolateOp());
      }
      else if (!tbox::Pointer<pdat::SideVariable<float> >(var).isNull()) {
         time_op = tbox::Pointer<hier::TimeInterpolateOperator>(
            new pdat::SideFloatLinearTimeInterpolateOp());
      }
      else {
         TBOX_ERROR("Unsupported variable centering for time interpolate operator");
      }
   }
   else {
      TBOX_ERROR("Unknown time interpolate operator name");
   }
   addTimeInterpolateOperator(time_op);
   return time_op;
}

}
}
#endif
