/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex data arrays.
 *
 ************************************************************************/

namespace SAMRAI {
namespace math {

SAMRAI_INLINE_KEYWORD
double
ArrayDataNormOpsComplex::L2NormWithControlVolume(
   const pdat::ArrayData<dcomplex>& data,
   const pdat::ArrayData<double>& cvol,
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(data, cvol, box);

   return sqrt(real(dotWithControlVolume(data, data, cvol, box)));
}

SAMRAI_INLINE_KEYWORD
double
ArrayDataNormOpsComplex::L2Norm(
   const pdat::ArrayData<dcomplex>& data,
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(data, box);

   return sqrt(real(dot(data, data, box)));
}

}
}
