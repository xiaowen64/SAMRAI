/*
 * File:	$URL$
 * Package:	SAMRAI patch data
 * Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:	$LastChangedRevision: 2244 $
 * Modified:	$LastChangedDate: 2008-07-02 11:10:50 -0700 (Wed, 02 Jul 2008) $
 * Description:	Light weight array class.
 */

#ifndef included_ArrayDataAccess_h
#define included_ArrayDataAccess_h

#ifndef included_tbox_ArrayData
#include "ArrayData.h"
#endif

#ifndef included_tbox_MDA_Access
#include "MDA_Access.h"
#endif

namespace SAMRAI {
namespace pdat {


/*!
 * @brief Utility for wrapping data from ArrayData class in an
 * MDA_Access or MDA_AccessConst object.
 *
 * This, and other classes in the folder array_access are meant
 * only for certain SAMRAI development efforts.  Please do not
 * use in user code.
 */
class ArrayDataAccess {

public:

   /*!
    * @brief Create an MDA_Access version of an ArrayData object.
    */
   template<int DIM, class TYPE>
   static MDA_Access<TYPE,DIM,MDA_OrderColMajor<DIM> >
   access( ArrayData<DIM,TYPE> &array_data, int depth=0 ) {
      return MDA_Access<TYPE,DIM,MDA_OrderColMajor<DIM> >(
         array_data.getPointer(depth),
         (const int*)array_data.getBox().lower(),
         (const int*)array_data.getBox().upper() );
   }

   /*!
    * @brief Create an MDA_AccessConst version of a const ArrayData object.
    */
   template<int DIM, class TYPE>
   static MDA_AccessConst<TYPE,DIM,MDA_OrderColMajor<DIM> >
   access( const ArrayData<DIM,TYPE> &array_data, int depth=0 ) {
      return MDA_AccessConst<TYPE,DIM,MDA_OrderColMajor<DIM> >(
         array_data.getPointer(depth),
         (const int*)array_data.getBox().lower(),
         (const int*)array_data.getBox().upper() );
   }

};

}
}

#endif
