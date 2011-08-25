/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating edge data objects
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int EdgeDataFactory<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
hier::MultiblockDataTranslator *
EdgeDataFactory<TYPE>::getMultiblockDataTranslator()
{
   if (d_mb_trans == NULL) {
      d_mb_trans = new MultiblockEdgeDataTranslator<TYPE>();
   }
   return d_mb_trans;
}

}
}
