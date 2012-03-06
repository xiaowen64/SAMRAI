/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex cell-centered patch data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace math {

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsComplex::printData(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsComplex::copyData(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& dst,
   const boost::shared_ptr<pdat::CellData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   dst->getArrayData().copy(src->getArrayData(), box);
}

SAMRAI_INLINE_KEYWORD
void
PatchCellDataOpsComplex::setToScalar(
   const boost::shared_ptr<pdat::CellData<dcomplex> >& dst,
   const dcomplex& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
