/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating node data objects 
 *
 ************************************************************************/

namespace SAMRAI {
namespace pdat {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int NodeDataFactory<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
hier::MultiblockDataTranslator *
NodeDataFactory<TYPE>::getMultiblockDataTranslator()
{
   if (d_mb_trans == NULL) {
      d_mb_trans = new MultiblockNodeDataTranslator<TYPE>();
   }
   return d_mb_trans;
}

}
}
