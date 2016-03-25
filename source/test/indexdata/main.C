#include <iostream>
#include <cassert>

#include "IndexVariable.h"
#include "IndexVariable.C"
#include "tbox/Array.h"
#include "tbox/Array.C"
#include "tbox/List.h"
#include "tbox/List.C"
#include "tbox/Pointer.h"
#include "tbox/Pointer.C"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/SAMRAIManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Timer.h"
#include "CellData.C"
#include "CellGeometry.C"
#include "IndexData.h"
#include "IndexData.C"
#include "IndexDataFactory.h"
#include "IndexDataFactory.C"

using namespace SAMRAI;
using namespace hier;
using namespace pdat;
using namespace tbox;

using namespace std;

#define NN 10

class Item {

public:
   Item()
   {
   };

   ~Item()
   {
   };
   
   void copySourceItem(const hier::Index<2>& idx,
                       const hier::IntVector<2>& src_offset,
                       const Item& src_item)
   {
      for (int n = 0; n < NN; n++) {
	 x[n] = src_item.x[n];
      }
   };
   
   size_t getDataStreamSize() { return 0; };
   void packStream(AbstractStream& stream) {};
   void unpackStream(AbstractStream& stream,
                     const hier::IntVector<2> offset) {};
   void putToDatabase(tbox::Pointer<tbox::Database> dbase) {};
   void getFromDatabase(tbox::Pointer<tbox::Database> dbase) {};

   double x[NN];
};

template class pdat::IndexData<2, Item, pdat::CellGeometry<2> >;
template class pdat::IndexDataFactory<2, Item, pdat::CellGeometry<2> >;
template class pdat::IndexDataNode<2, Item, pdat::CellGeometry<2> >;
template class pdat::IndexIterator<2, Item, pdat::CellGeometry<2> >;
template class pdat::IndexVariable<2, Item, pdat::CellGeometry<2> >;
template class tbox::Array< IndexDataNode<2, Item, pdat::CellGeometry<2> > >;

template class tbox::List< pdat::IndexDataNode<2, Item, pdat::CellGeometry<2> > >;
template class tbox::ListIterator< pdat::IndexDataNode<2, Item, pdat::CellGeometry<2> > >;
template class tbox::ListNode< pdat::IndexDataNode<2, Item, pdat::CellGeometry<2> > >;
template class tbox::Pointer< pdat::IndexData<2, Item, pdat::CellGeometry<2> > >;
template class tbox::Pointer< pdat::IndexVariable<2, Item, pdat::CellGeometry<2> > >;
template class tbox::Pointer< IndexDataFactory<2, Item, pdat::CellGeometry<2> > >;

int main(int argc, char* argv[])
{
   SAMRAI_MPI::init(&argc,&argv);
   SAMRAIManager::startup();

   
   Index<2> lo = Index<2>(0);
   Index<2> hi = Index<2>(100);
   Box<2> box(lo,hi);
   
   srand(1);

   /******************************************************************************
    * InedxData interface tests.
    ******************************************************************************/
   {
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      Item* item = new Item;
      Index<2> idx(0,0);
      idx_data.addItemPointer(idx, item);

      // isElement()
      assert( idx_data.isElement(idx) );
      Index<2> idx2(1,0);
      assert( !idx_data.isElement(idx2) );
      Index<2> idx3(0,1);
      assert( !idx_data.isElement(idx3) );

      // addItem()/getItem()
      assert( idx_data.getItem(idx) == item );

      assert( idx_data.getNumberOfItems() == 1);
      
      // removeItem()
      idx_data.removeItem(idx);
      assert( !idx_data.isElement(idx) );

      assert( idx_data.getNumberOfItems() == 0);
   }

   {
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      Item* item = new Item;
      Index<2> idx(0,0);
      idx_data.addItem(idx, *item);
      delete item;

      // isElement()
      assert( idx_data.isElement(idx) );
      Index<2> idx2(1,0);
      assert( !idx_data.isElement(idx2) );
      Index<2> idx3(0,1);
      assert( !idx_data.isElement(idx3) );

      // addItem()/getItem()
      assert( idx_data.getNumberOfItems() == 1);
      
      // removeItem()
      idx_data.removeItem(idx);
      assert( !idx_data.isElement(idx) );

      assert( idx_data.getNumberOfItems() == 0);
   }

   {
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      Item* item = new Item;
      Index<2> idx(0,0);
      idx_data.replaceAddItem(idx, *item);
      delete item;

      // isElement()
      assert( idx_data.isElement(idx) );
      Index<2> idx2(1,0);
      assert( !idx_data.isElement(idx2) );
      Index<2> idx3(0,1);
      assert( !idx_data.isElement(idx3) );

      // addItem()/getItem()
      assert( idx_data.getNumberOfItems() == 1);

      item = new Item;
      idx_data.replaceAddItem(idx, *item);
      delete item;

      assert( idx_data.getNumberOfItems() == 1);
      
      // removeItem()
      idx_data.removeItem(idx);
      assert( !idx_data.isElement(idx) );
   }


   {
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      // getNumberItems()

      Index<2> idx1(0,0);
      idx_data.addItemPointer(idx1, new Item);

      Index<2> idx2(1,0);
      idx_data.addItemPointer(idx2, new Item);

      assert( idx_data.getNumberOfItems() == 2 );

      // remove 1
      idx_data.removeItem(idx1);
      assert( idx_data.getNumberOfItems() == 1 );

      // replace 1 at same index, no change
      idx_data.addItemPointer(idx2, new Item);
      assert( idx_data.getNumberOfItems() == 1 );
   }

   {
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      // removeInsideBox()
      Index<2> lo(2,2);
      Index<2> hi(3,5);

      Box<2> box1(lo,hi);
      for (Box<2>::Iterator bi(box1); bi; bi++) {

         Index<2> idx = bi();

         idx_data.addItemPointer(idx, new Item);

      }

      assert( idx_data.getNumberOfItems() == box1.size() );

      idx_data.removeInsideBox(box1);

      assert (idx_data.getNumberOfItems() == 0 );
   }
   {
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      // removeAllItems()
      Index<2> lo(0,0);
      Index<2> hi(1,1);

      Box<2> box1(lo,hi);
      for (Box<2>::Iterator bi(box1); bi; bi++) {

         Index<2> idx = bi();

         idx_data.addItemPointer(idx, new Item);

      }

      assert( idx_data.getNumberOfItems() == box1.size() );

      idx_data.removeAllItems();

      assert (idx_data.getNumberOfItems() == 0 );
   }

   {
      // copy() where src and dst are same box
      
      IndexData<2, Item, pdat::CellGeometry<2> > src(box, 0);
      IndexData<2, Item, pdat::CellGeometry<2> > dst(box, 0);


      Index<2> lo(0,0);
      Index<2> hi(1,1);

      Box<2> box1(lo,hi);
      for (Box<2>::Iterator bi(box1); bi; bi++) {
         src.addItemPointer(bi(), new Item);
      }

      assert( src.getNumberOfItems() == box1.size() );
      assert( dst.getNumberOfItems() == 0 );
      
      dst.copy(src);

      assert( dst.getNumberOfItems() == src.getNumberOfItems() );
   }

   {

      // copy() where src and dst partially overlap, and only
      // some of src's items are contained in overlap.

      Index<2> lo_src(0,0);
      Index<2> hi_src(2,2);
      Box<2> box_src(lo_src,hi_src);
      IndexData<2, Item, pdat::CellGeometry<2> > src(box_src, 0);

      // Two of these three items should end up in dst
      Index<2> idx_item1(0,0);
      src.addItemPointer(idx_item1, new Item);

      Index<2> idx_item2(1,1);
      src.addItemPointer(idx_item2, new Item);

      Index<2> idx_item3(2,2);
      src.addItemPointer(idx_item3, new Item);

      
      Index<2> lo_dst(1,1);
      Index<2> hi_dst(3,3);
      Box<2> box_dst(lo_dst,hi_dst);

      IndexData<2, Item, pdat::CellGeometry<2> > dst(box_dst, 0);

      assert( src.getNumberOfItems() == 3 );
      assert( dst.getNumberOfItems() == 0 );
      
      dst.copy(src);

      assert( dst.getNumberOfItems() == 2 );
   }

   {
      // copy() with overlap argument

      // src
      // x . . . . . x "
      // .   .   . 3 . " // 2
      // . . . . . . . "
      // .   . 2 .   . " // 1
      // . . . . . . . "
      // . 1 .   .   . " // 0
      // x . . . . . x "

      // dst orig
      // . . x . . . x "
      // .   . 4 .   . " // 2
      // . . . . . . . "
      // .   .   .   . " // 1
      // . . x . . . x "
      // .   .   .   . " // 0
      // . . . . . . . "
      
      // dst expected
      // . . x . . . x "
      // .   .   . 3 . "
      // . . . . . . . "
      // .   . 2 .   . "
      // . . x . . . x "
      // .   .   .   . "
      // . . . . . . . "
      
      Index<2> lo_src(0,0);
      Index<2> hi_src(2,2);
      Box<2> box_src(lo_src,hi_src);
      IndexData<2, Item, pdat::CellGeometry<2> > src(box_src, 0);

      // Two of these three items should end up in dst
      Index<2> idx_item1(0,0);
      src.addItemPointer(idx_item1, new Item);

      Index<2> idx_item2(1,1);
      src.addItemPointer(idx_item2, new Item);

      Index<2> idx_item3(2,2);
      src.addItemPointer(idx_item3, new Item);
      
      Index<2> lo_dst(1,1);
      Index<2> hi_dst(2,2);
      Box<2> box_dst(lo_dst,hi_dst);

      IndexData<2, Item, pdat::CellGeometry<2> > dst(box_dst, 0);

      // This item should be removed
      Index<2> idx_item4(1,2);
      dst.addItemPointer(idx_item4, new Item);

      assert( src.getNumberOfItems() == 3 );
      assert( dst.getNumberOfItems() == 1 );

      IntVector<2> src_offset(0);
      BoxList<2> boxes(box_src);
      boxes.addItem(box_dst);
      BoxList<2> intersection = box_src * box_dst;
      CellOverlap<2> overlap(intersection, src_offset);

      BoxList<2> dst_boxlist(overlap.getDestinationBoxList());

      dst.copy(src, overlap);

      assert( dst.getNumberOfItems() == 2 );
   }

   {
      // copy(): Same as test7 using copy2 which reverses src and dst

      Index<2> lo_src(0,0);
      Index<2> hi_src(2,2);
      Box<2> box_src(lo_src,hi_src);
      IndexData<2, Item, pdat::CellGeometry<2> > src(box_src, 0);

      // Two of these three items should end up in dst
      Index<2> idx_item1(0,0);
      src.addItemPointer(idx_item1, new Item);

      Index<2> idx_item2(1,1);
      src.addItemPointer(idx_item2, new Item);

      Index<2> idx_item3(2,2);
      src.addItemPointer(idx_item3, new Item);

      
      Index<2> lo_dst(1,1);
      Index<2> hi_dst(3,3);
      Box<2> box_dst(lo_dst,hi_dst);

      IndexData<2, Item, pdat::CellGeometry<2> > dst(box_dst, 0);

      assert( src.getNumberOfItems() == 3 );
      assert( dst.getNumberOfItems() == 0 );

      src.copy2(dst);

      assert( dst.getNumberOfItems() == 2 );
   }

   {
      Index<2> lo(0,0);
      Index<2> hi(2,2);
      Box<2> box(lo,hi);
      IndexData<2, Item, pdat::CellGeometry<2> > data(box, 0);

      // Add three items
      Index<2> idx_item1(0,0);
      data.addItemPointer(idx_item1, new Item);

      Index<2> idx_item2(0,1);
      data.addItemPointer(idx_item2, new Item);

      Index<2> idx_item3(2,1);
      data.addItemPointer(idx_item3, new Item);

      int count = 0;
      for (IndexIterator<2, Item, pdat::CellGeometry<2> > it(data); it; it++) {
         count++;
      }
      assert(3 == count);
   }
   
   int size=100;
   {
      tbox::Pointer<tbox::Timer> timer;

      timer =   tbox::TimerManager::getManager() ->
	 getTimer("IndexDataAppendItemSequential", true);
      
      tbox::pout << "Begin Timing" << endl;
      
      Index<2> lo = Index<2>(0);
      Index<2> hi = Index<2>(size);
      Box<2> box(lo,hi);
   
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      timer -> start();

      for(int i = 0; i < size; ++i) {
	 for(int j = 0; j < size; ++j) {
	    Index<2> idx(i,j);
	    
	    Item new_item;
	    idx_data.appendItem(idx, new_item);
	 }
      }

      int numberOfItems = idx_data.getNumberOfItems();
      timer -> stop();
      
      tbox::pout << numberOfItems << endl;

      tbox::pout.precision(16);

      tbox::pout << "IndexData appendItem Sequential insert time : " << timer -> getTotalWallclockTime() << endl;

      tbox::pout << "End Timing" << endl;
   }

   {
      tbox::Pointer<tbox::Timer> timer;

      timer =   tbox::TimerManager::getManager() ->
	 getTimer("IndexDataAppendItemPointerSequential", true);

      tbox::pout << "Begin Timing" << endl;

      Index<2> lo = Index<2>(0);
      Index<2> hi = Index<2>(size);
      Box<2> box(lo,hi);
   
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      timer -> start();

      for(int i = 0; i < size; ++i) {
	 for(int j = 0; j < size; ++j) {
	    Index<2> idx(i,j);
	    
	    Item *new_item = new Item();
	    idx_data.appendItemPointer(idx, new_item);
	 }
      }

      timer -> stop();

      tbox::pout.precision(16);
      
      tbox::pout << "IndexData appendItemPointer sequential insert time : " << timer -> getTotalWallclockTime() << endl;

      tbox::pout << "End Timing" << endl;
   }


   size = 100;
   int num_inserts = 1e5;

   {
      tbox::Pointer<tbox::Timer> timer;

      timer =   tbox::TimerManager::getManager() ->
	 getTimer("IndexDataAppendItemRandom", true);
      
      tbox::pout << "Begin Timing" << endl;
      
      Index<2> lo = Index<2>(0);
      Index<2> hi = Index<2>(size);
      Box<2> box(lo,hi);
   
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      timer -> start();

      for (int n = 0; n < num_inserts; n++) {
	 int i = rand() % size;
	 int j = rand() % size;
	 Index<2> idx(i,j);
	    
	 Item new_item;
	 idx_data.appendItem(idx, new_item);
      }

      int numberOfItems = idx_data.getNumberOfItems();
      timer -> stop();
      
      tbox::pout << numberOfItems << endl;

      tbox::pout.precision(16);

      tbox::pout << "IndexData appendItem random insert time : " << timer -> getTotalWallclockTime() << endl;

      tbox::pout << "End Timing" << endl;
   }

   {
      tbox::Pointer<tbox::Timer> timer;

      timer =   tbox::TimerManager::getManager() ->
	 getTimer("IndexDataAppendItemPointerRandom", true);

      tbox::pout << "Begin Timing" << endl;

      Index<2> lo = Index<2>(0);
      Index<2> hi = Index<2>(size);
      Box<2> box(lo,hi);
   
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      timer -> start();

      for (int n = 0; n < num_inserts; n++) {
	 int i = rand() % size;
	 int j = rand() % size;
	 Index<2> idx(i,j);
	    
	 Item *new_item = new Item();
	 idx_data.appendItemPointer(idx, new_item);
      }
      
      timer -> stop();

      tbox::pout.precision(16);
      
      tbox::pout << "IndexData appendItemPointer random insert time : " << timer -> getTotalWallclockTime() << endl;

      tbox::pout << "End Timing" << endl;
   }

   size = 100;

   {
      tbox::Pointer<tbox::Timer> timer;

      timer =   tbox::TimerManager::getManager() ->
	 getTimer("IndexDataReplace", true);

      tbox::pout << "Begin Timing" << endl;

      Index<2> lo = Index<2>(0);
      Index<2> hi = Index<2>(size);
      Box<2> box(lo,hi);
   
      IndexData<2, Item, pdat::CellGeometry<2> > idx_data(box, 0);

      timer -> start();

      for (int n = 0; n < num_inserts; n++) {
	 int i = rand() % size;
	 int j = rand() % size;
	 Index<2> idx(i,j);
	    
	 Item *new_item = new Item();
	 idx_data.replaceAddItemPointer(idx, new_item);
      }
      
      timer -> stop();

      tbox::pout.precision(16);
      
      tbox::pout << "IndexData replaceAddItemPointer random insert time : " << timer -> getTotalWallclockTime() << endl;

      tbox::pout << "End Timing" << endl;
   }

   tbox::pout << "PASSED" << endl;

   SAMRAIManager::shutdown();
   SAMRAI_MPI::finalize();

   exit(0);
}

