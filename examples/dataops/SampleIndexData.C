//
// File:        SampleIndexData.C
// Package:     SAMRAI application
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 480 $
// Modified:    $Date: 2005-07-14 13:46:02 -0700 (Thu, 14 Jul 2005) $
// Description: SampleIndexData example demonstrating IndexData type.
//

#include "SampleIndexData.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#include "tbox/Utilities.h"


/*
*************************************************************************
*                                                                       *
* Default Constructor                                                   *
*                                                                       *
*************************************************************************
*/

SampleIndexData::SampleIndexData()
{
   d_dummy_int = 0;
}

/*
*************************************************************************
*                                                                       *
* Constructor providing cell index.                                     *
*                                                                       *
*************************************************************************
*/

SampleIndexData::SampleIndexData(const pdat::CellIndex<NDIM>& ic)
{
   d_index = ic;
   d_dummy_int = 0;
}


/*
*************************************************************************
*                                                                       *
* Copy Constructor                                                     *
*                                                                       *
*************************************************************************
*/

SampleIndexData::SampleIndexData(const SampleIndexData& data) :
   d_index(data.d_index),
   d_dummy_int(data.d_dummy_int)
{

}

/*
*************************************************************************
*                                                                       *
* Assignment operator                                                   *
*                                                                       *
*************************************************************************
*/

SampleIndexData& SampleIndexData::operator=(const SampleIndexData& data)
{
   d_index = data.d_index;
   d_dummy_int = data.d_dummy_int;
   return(*this);
}

/*
*************************************************************************
*                                                                       *
* Destructor
*                                                                       *
*************************************************************************
*/

SampleIndexData::~SampleIndexData()
{
}

/*
*************************************************************************
*  
* Set dummy int data
*                                                                       *
*************************************************************************
*/
void SampleIndexData::setInt(const int dummy) 
{
   d_dummy_int = dummy;
}

/*
*************************************************************************
*  
*  Return dummy int data
*                                                                       *
*************************************************************************
*/
int SampleIndexData::getInt() const 
{
   return(d_dummy_int);
}

/*
*************************************************************************
*  
*  Return index
*                                              *
*************************************************************************
*/
pdat::CellIndex<NDIM> SampleIndexData::getIndex() const 
{
   return(d_index);
}

/*
*************************************************************************
*  
*  Print contents
*                                              *
*************************************************************************
*/
void SampleIndexData::printClassData(ostream& os) const 
{
   
}

/*
*************************************************************************
*   
* The copySourceItem() method allows SampleIndexData to be a templated
* data type for IndexData - i.e. IndexData<SampleIndexData>.  
*                                                                       *
*************************************************************************
*/
void SampleIndexData::copySourceItem(hier::Index<NDIM>& index,
                                const hier::IntVector<NDIM>& src_offset,
                                SampleIndexData& src_item)
{

   /*
    * Copy src_item data into *this.  Note that we don't do
    * anything with the src_offset.  This is because we have
    * access to the index already.
    */
   d_index           = (pdat::CellIndex<NDIM>) index;
   d_dummy_int    = src_item.d_dummy_int;
   

}

/*
*************************************************************************
*                                                                       *
* The getDataStreamSize(), packStream(), and unpackStream() methods
* are required to template SampleIndexData as IndexData type - i.e.
* IndexData<SampleIndexData>.  They are used to communicate SampleIndexData,
* specifying how many bytes will be packed during the "packStream()" 
* method.
*                                                                       *
*************************************************************************
*/

size_t SampleIndexData::getDataStreamSize()
{
   /*
    * #bytes = 
    *   d_index           (int[NDIM]) +
    *   d_dummy_int       (int) 
    */
   size_t bytes = (NDIM+1)*tbox::AbstractStream::sizeofInt();

   return(bytes);
}

void SampleIndexData::packStream(tbox::AbstractStream& stream)
{
   int counter = 0;
   int ibuffer[NDIM+1];
   for (int i=0; i<NDIM; i++) {
      ibuffer[i] = d_index(i);
      counter++;
   }
   ibuffer[counter] = d_dummy_int;   
   stream.pack(ibuffer, NDIM+1);

}

void SampleIndexData::unpackStream(tbox::AbstractStream& stream,
                                const hier::IntVector<NDIM>& offset)
{
   int counter = 0;
   int ibuffer[NDIM+1];
   stream.unpack(ibuffer, NDIM);
   pdat::CellIndex<NDIM> index;   
   for (int i=0; i<NDIM; i++) {
      index[i] = ibuffer[i];
      counter++;
   }
   d_index = index + offset;

   d_dummy_int = ibuffer[counter];   

}



/*
*************************************************************************
*                                                                       *
* The putToDatabase() and getFromDatabase() methods
* are required to template SampleIndexData as IndexData type - i.e.
* IndexData<SampleIndexData>.  They are used to write/read SampleIndexData,
* data to/from the restart database. 
*                                                                       *
*************************************************************************
*/

void SampleIndexData::putToDatabase(
   tbox::Pointer<tbox::Database>& database)
{

   int counter = 0;
   int ibuffer[NDIM+1];
   for (int i=0; i<NDIM; i++) {
      ibuffer[i] = d_index(i);
      counter++;
   }
   ibuffer[counter] = d_dummy_int;
   
   database->putIntegerArray("ibuffer", ibuffer, NDIM+1);

}


void SampleIndexData::getFromDatabase(
   tbox::Pointer<tbox::Database>& database)
{
   int ibuffer[NDIM+1];
   database->getIntegerArray("ibuffer", ibuffer, NDIM+1);
   pdat::CellIndex<NDIM> index;
   for (int i=0; i<NDIM; i++) {
      index(i) = ibuffer[i];
   }
   d_index = index;
   d_dummy_int = ibuffer[NDIM];
   
}

/*
*****************************************************************
* 
*  Templates used for SampleIndexData
*
*****************************************************************
*/

#include "SampleIndexData.h"
#include "tbox/Array.C"
#include "IndexData.C"
#include "IndexDataFactory.C"
#include "IndexVariable.C"
#include "tbox/Pointer.C"

namespace SAMRAI {

template class pdat::IndexData<NDIM, SampleIndexData >;
template class pdat::IndexDataFactory<NDIM,SampleIndexData>;
template class pdat::IndexDataNode<NDIM,SampleIndexData>;
template class pdat::IndexIterator<NDIM, SampleIndexData >;
template class pdat::IndexVariable<NDIM, SampleIndexData >;
template class tbox::Array< SampleIndexData >;
template class tbox::Array< pdat::IndexDataNode<NDIM, SampleIndexData > >;
template class tbox::Pointer < pdat::IndexData<NDIM, SampleIndexData > >;
template class tbox::Pointer < pdat::IndexVariable<NDIM, SampleIndexData > >;
template class tbox::Pointer < pdat::IndexDataFactory<NDIM, SampleIndexData > >;

}
