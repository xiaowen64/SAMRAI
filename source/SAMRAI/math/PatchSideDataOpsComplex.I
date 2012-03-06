/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex side-centered patch data.
 *
 ************************************************************************/

namespace SAMRAI {
namespace math {

SAMRAI_INLINE_KEYWORD
void
PatchSideDataOpsComplex::printData(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& data,
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
PatchSideDataOpsComplex::copyData(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& dst,
   const boost::shared_ptr<pdat::SideData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         dst->getArrayData(d).copy(src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}

SAMRAI_INLINE_KEYWORD
void
PatchSideDataOpsComplex::setToScalar(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& dst,
   const dcomplex& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
