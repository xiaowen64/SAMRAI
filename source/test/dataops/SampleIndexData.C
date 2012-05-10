/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   SampleIndexData example demonstrating IndexData type.
 *
 ************************************************************************/

#include "SampleIndexData.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Dimension.h"

#include <iostream>

/*
 *************************************************************************
 *
 * Constructor providing cell index.
 *
 *************************************************************************
 */

using namespace SAMRAI;

SampleIndexData::SampleIndexData():
   d_dummy_int(0)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

SampleIndexData::~SampleIndexData()
{
}

/*
 *************************************************************************
 *
 * Set dummy int data
 *
 *************************************************************************
 */
void SampleIndexData::setInt(
   const int dummy)
{
   d_dummy_int = dummy;
}

/*
 *************************************************************************
 *
 *  Return dummy int data
 *
 *************************************************************************
 */
int SampleIndexData::getInt() const
{
   return d_dummy_int;
}

/*
 *************************************************************************
 *
 * The copySourceItem() method allows SampleIndexData to be a templated
 * data type for IndexData - i.e. IndexData<SampleIndexData>.
 *
 *************************************************************************
 */
void SampleIndexData::copySourceItem(
   const hier::Index& index,
   const hier::IntVector& src_offset,
   const SampleIndexData& src_item)
{
   NULL_USE(index);
   NULL_USE(src_offset);
   d_dummy_int = src_item.d_dummy_int;
}

/*
 *************************************************************************
 *
 * The getDataStreamSize(), packStream(), and unpackStream() methods
 * are required to template SampleIndexData as IndexData type - i.e.
 * IndexData<SampleIndexData>.  They are used to communicate SampleIndexData,
 * specifying how many bytes will be packed during the "packStream()"
 * method.
 *
 *************************************************************************
 */

size_t SampleIndexData::getDataStreamSize()
{
   return 0;
}

void SampleIndexData::packStream(
   tbox::MessageStream& stream)
{
   NULL_USE(stream);
}

void SampleIndexData::unpackStream(
   tbox::MessageStream& stream,
   const hier::IntVector& offset)
{
   NULL_USE(stream);
   NULL_USE(offset);
}

/*
 *************************************************************************
 *
 * The putToDatabase() and getFromDatabase() methods
 * are required to template SampleIndexData as IndexData type - i.e.
 * IndexData<SampleIndexData>.  They are used to write/read SampleIndexData,
 * data to/from the restart database.
 *
 *************************************************************************
 */

void SampleIndexData::putUnregisteredToDatabase(
   boost::shared_ptr<tbox::Database>& database) const
{
   NULL_USE(database);
}

void SampleIndexData::getFromDatabase(
   boost::shared_ptr<tbox::Database>& database)
{
   NULL_USE(database);
}

/*
 *****************************************************************
 *
 *  Templates used for SampleIndexData
 *
 *****************************************************************
 */

//#include "SampleIndexData.h"
//#include "SAMRAI/tbox/Array.C"
//#include "SAMRAI/pdat/IndexData.C"
//#include "SAMRAI/pdat/IndexDataFactory.C"
//#include "SAMRAI/pdat/IndexVariable.C"
//#include "SAMRAI/pdat/CellGeometry.h"
//
//namespace SAMRAI {
//
//template class pdat::SparseData<SampleIndexData, pdat::CellGeometry>;
//template class pdat::SparseDataFactory<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexData<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexDataFactory<SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexDataNode<NDIM, SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexIterator<NDIM, SampleIndexData, pdat::CellGeometry>;
//template class pdat::IndexVariable<SampleIndexData, pdat::CellGeometry>;
//template class tbox::Array<SampleIndexData>;
//template class tbox::Array<pdat::IndexDataNode<NDIM, SampleIndexData,
//                                               pdat::CellGeometry> >;
//template class boost::shared_ptr<pdat::IndexData<NDIM, SampleIndexData,
//                                                 pdat::CellGeometry> >;
//template class boost::shared_ptr<pdat::IndexVariable<SampleIndexData,
//                                                     pdat::CellGeometry> >;
//template class boost::shared_ptr<pdat::IndexDataFactory<SampleIndexData,
//                                                        pdat::CellGeometry> >;
//
//}
