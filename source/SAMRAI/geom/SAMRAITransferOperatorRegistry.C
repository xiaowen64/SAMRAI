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

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace geom {

/*
 *************************************************************************
 *
 * Constructor and destructor for TransferOperatorRegistry objects.
 *
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

boost::shared_ptr<hier::CoarsenOperator>
SAMRAITransferOperatorRegistry::buildCoarsenOperator(
   const boost::shared_ptr<hier::Variable>& var,
   const std::string& op_name)
{
   boost::shared_ptr<hier::CoarsenOperator> coarsen_op;
   if (op_name == "CONSERVATIVE_COARSEN") {
      if (boost::shared_ptr<pdat::CellVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianCellComplexWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianCellDoubleWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianCellFloatWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianEdgeComplexWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianEdgeDoubleWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianEdgeFloatWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianFaceComplexWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianFaceDoubleWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianFaceFloatWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianOuterfaceComplexWeightedAverage>(
               getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianOuterfaceDoubleWeightedAverage>(
               getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianOuterfaceFloatWeightedAverage>(
               getDim());
      } else if (boost::shared_ptr<pdat::OutersideVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianOutersideDoubleWeightedAverage>(
               getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianSideComplexWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianSideDoubleWeightedAverage>(getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<CartesianSideFloatWeightedAverage>(getDim());
      } else {
         TBOX_ERROR("Unsupported variable centering for coarsen operator");
      }
   } else if (op_name == "SKELETON_COARSEN") {
      coarsen_op = boost::make_shared<SkeletonCoarsen>(getDim());
   } else if (op_name == "CONSTANT_COARSEN") {
      if (boost::shared_ptr<pdat::NodeVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op = boost::make_shared<pdat::NodeComplexInjection>(getDim());
      } else if (boost::shared_ptr<pdat::NodeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op = boost::make_shared<pdat::NodeDoubleInjection>(getDim());
      } else if (boost::shared_ptr<pdat::NodeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op = boost::make_shared<pdat::NodeFloatInjection>(getDim());
      } else if (boost::shared_ptr<pdat::NodeVariable<int> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op = boost::make_shared<pdat::NodeIntegerInjection>(getDim());
      } else if (boost::shared_ptr<pdat::OuternodeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         coarsen_op =
            boost::make_shared<pdat::OuternodeDoubleConstantCoarsen>(getDim());
      } else {
         TBOX_ERROR("Unsupported variable centering for coarsen operator");
      }
   } else {
      TBOX_ERROR("Unknown coarsen operator name");
   }
   addCoarsenOperator(coarsen_op);
   return coarsen_op;
}

boost::shared_ptr<hier::RefineOperator>
SAMRAITransferOperatorRegistry::buildRefineOperator(
   const boost::shared_ptr<hier::Variable>& var,
   const std::string& op_name)
{
   boost::shared_ptr<hier::RefineOperator> refine_op;
   if (op_name == "CONSERVATIVE_LINEAR_REFINE") {
      if (boost::shared_ptr<pdat::CellVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianCellComplexConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianCellDoubleConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianCellFloatConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianEdgeDoubleConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianEdgeFloatConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianFaceDoubleConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianFaceFloatConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianSideDoubleConservativeLinearRefine>(
               getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianSideFloatConservativeLinearRefine>(
               getDim());
      } else {
         TBOX_ERROR("Unsupported variable centering for refine operator");
      }
   } else if (op_name == "LINEAR_REFINE") {
      if (boost::shared_ptr<pdat::CellVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianCellComplexLinearRefine>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianCellDoubleLinearRefine>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianCellFloatLinearRefine>(getDim());
      } else if (boost::shared_ptr<pdat::NodeVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianNodeComplexLinearRefine>(getDim());
      } else if (boost::shared_ptr<pdat::NodeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianNodeDoubleLinearRefine>(getDim());
      } else if (boost::shared_ptr<pdat::NodeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<CartesianNodeFloatLinearRefine>(getDim());
      } else {
         TBOX_ERROR("Unsupported variable centering for refine operator");
      }
   } else if (op_name == "SKELETON_REFINE") {
      refine_op = boost::make_shared<SkeletonRefine>(getDim());
   } else if (op_name == "CONSTANT_REFINE") {
      if (boost::shared_ptr<pdat::CellVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::CellComplexConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::CellDoubleConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::CellFloatConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::CellVariable<int> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::CellIntegerConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::EdgeComplexConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::EdgeDoubleConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::EdgeFloatConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::EdgeVariable<int> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::EdgeIntegerConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::FaceComplexConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::FaceDoubleConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::FaceFloatConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::FaceVariable<int> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::FaceIntegerConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::OuterfaceComplexConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::OuterfaceDoubleConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::OuterfaceFloatConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<int> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::OuterfaceIntegerConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::SideComplexConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::SideDoubleConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::SideFloatConstantRefine>(getDim());
      } else if (boost::shared_ptr<pdat::SideVariable<int> >(var, boost::detail::dynamic_cast_tag())) {
         refine_op =
            boost::make_shared<pdat::SideIntegerConstantRefine>(getDim());
      } else {
         TBOX_ERROR("Unsupported variable centering for refine operator");
      }
   } else {
      TBOX_ERROR("Unknown coarsen operator name");
   }
   addRefineOperator(refine_op);
   return refine_op;
}

boost::shared_ptr<hier::TimeInterpolateOperator>
SAMRAITransferOperatorRegistry::buildTimeInterpolateOperator(
   const boost::shared_ptr<hier::Variable>& var,
   const std::string& op_name)
{
   boost::shared_ptr<hier::TimeInterpolateOperator> time_op;
   if (op_name == "STD_LINEAR_TIME_INTERPOLATE") {
      if (boost::shared_ptr<pdat::CellVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::CellComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::CellVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::CellDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::CellVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::CellFloatLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::EdgeVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::EdgeComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::EdgeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::EdgeDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::EdgeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::EdgeFloatLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::FaceVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::FaceComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::FaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::FaceDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::FaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::FaceFloatLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::NodeVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::NodeComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::NodeVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::NodeDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::NodeVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::NodeFloatLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::OuterfaceComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::OuterfaceDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::OuterfaceVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::OuterfaceFloatLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::OutersideVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::OutersideComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::OutersideVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::OutersideDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::OutersideVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::OutersideFloatLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::SideVariable<dcomplex> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::SideComplexLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::SideVariable<double> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::SideDoubleLinearTimeInterpolateOp>();
      } else if (boost::shared_ptr<pdat::SideVariable<float> >(var, boost::detail::dynamic_cast_tag())) {
         time_op =
            boost::make_shared<pdat::SideFloatLinearTimeInterpolateOp>();
      } else {
         TBOX_ERROR("Unsupported variable centering for time interpolate operator");
      }
   } else {
      TBOX_ERROR("Unknown time interpolate operator name");
   }
   addTimeInterpolateOperator(time_op);
   return time_op;
}

}
}
#endif
