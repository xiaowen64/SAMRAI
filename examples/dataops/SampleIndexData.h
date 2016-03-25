//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/dataops/SampleIndexData.h $
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Boundary cell struct for embedded boundary implementations
//

#ifndef included_SampleIndexDataXD
#define included_SampleIndexDataXD

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "tbox/AbstractStream.h"
#include "CellIndex.h"
#include "tbox/Database.h"
#include "Index.h"
#include "IntVector.h"
#include "tbox/Pointer.h"
#include "tbox/IOStream.h"


/**
 * The SampleClass struct holds some dummy data and methods.  It's intent
 * is to indicate how a user could construct their own index data type.
 */

using namespace SAMRAI;

class SampleIndexData
{
public:
   /**
    * The default constructor creates an ``empty'' SampleIndexData.
    */
   SampleIndexData();

   /**
    * Copy constructor.
    */
   SampleIndexData(const SampleIndexData& data);

   /**
    * Constructor supplying cell index where data is defined.
    */
   SampleIndexData(const pdat::CellIndex<NDIM>& ic);

   /**
    * The assignment operator copies the data of the argument cell.
    */
   SampleIndexData& operator=(const SampleIndexData& cell);

   /**
    * The destructor for SampleIndexData.
    */
   ~SampleIndexData();

   /**
    * Sets a dummy integer in this class.
    */
   void setInt(const int dummy);

   /**
    * Returns a dummy integer in this class.
    */
   int getInt() const;

   /**
    * Returns the cell index where the index data is stored.
    */
   pdat::CellIndex<NDIM> getIndex() const;

   /**
    * Print class data representation when an unrecoverable run-time
    * exception is thrown. Or, when desired.
    */
   void printClassData(ostream& os) const;

   /**
    * The copySourceItem() method allows SampleIndexData to be a templated
    * data type for IndexData - i.e. IndexData<SampleIndexData>.  In
    * addition to this method, the other methods that must be defined are
    * getDataStreamSize(), packStream(), unpackStream() for communication,
    * putToDatabase(), getFromDatabase for restart.  These are described
    * below. 
    */
   void copySourceItem(hier::Index<NDIM>& index, 
                       const hier::IntVector<NDIM>& src_offset, 
                       SampleIndexData& src_item);

   /**
    * The following functions enable parallel communication with SampleIndexDatas.
    * They are used in SAMRAI communication infrastructure to 
    * specify the number of bytes of data stored in each SampleIndexData object,
    * and to pack and unpack the data to the specified stream.
    */
   size_t getDataStreamSize();
   void packStream(tbox::AbstractStream& stream);
   void unpackStream(tbox::AbstractStream& stream, const hier::IntVector<NDIM>& offset);
   

   /**
    * These functions are used to read/write SampleIndexData data to/from 
    * restart.
    */
   void getFromDatabase(tbox::Pointer<tbox::Database>& database);
   void putToDatabase(tbox::Pointer<tbox::Database>& database);
   
 
private:
   /*
    * Cell index where SampleIndexData is defined.
    */
   pdat::CellIndex<NDIM> d_index;
   
   /*
    * Dummy int data
    */
   int d_dummy_int;

   // ADD ANY OTHER DATA HERE
};
#endif 
