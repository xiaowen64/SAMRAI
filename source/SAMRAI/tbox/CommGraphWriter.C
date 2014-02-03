/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Collects and writes out data on communication graphs.
 *
 ************************************************************************/
#include "SAMRAI/tbox/CommGraphWriter.h"


#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif


#include "SAMRAI/tbox/MessageStream.h"

namespace SAMRAI {
namespace tbox {



/*
 ***********************************************************************
 ***********************************************************************
 */
CommGraphWriter::CommGraphWriter()
   : d_root_rank(0)
{
}



/*
 ***********************************************************************
 ***********************************************************************
 */
CommGraphWriter::~CommGraphWriter()
{
}



/*
 ***********************************************************************
 ***********************************************************************
 */
size_t CommGraphWriter::addRecord(
   const SAMRAI_MPI &mpi,
   size_t number_of_edge_types,
   size_t number_of_node_value_types )
{
   d_records.resize( 1 + d_records.size() );
   Record &record = d_records.back();
   record.d_mpi = mpi;
   record.d_edges.resize(number_of_edge_types);
   record.d_node_values.resize(number_of_node_value_types);
   return (d_records.size() - 1);
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::setEdgeInCurrentRecord(
   size_t edge_type_index,
   const std::string &edge_label,
   double edge_value,
   EdgeDirection edge_direction,
   int other_node )
{
   TBOX_ASSERT( edge_type_index < d_records.back().d_edges.size() );

   Edge &edge = d_records.back().d_edges[edge_type_index];

   edge.d_label = edge_label;
   edge.d_value = edge_value;
   edge.d_dir = edge_direction;
   edge.d_other_node = other_node;

   return;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::setNodeValueInCurrentRecord(
   size_t nodevalue_type_index,
   const std::string &nodevalue_label,
   double node_value )
{
   TBOX_ASSERT( nodevalue_type_index < d_records.back().d_node_values.size() );

   NodeValue &nodevalue = d_records.back().d_node_values[nodevalue_type_index];

   nodevalue.d_value = node_value;
   nodevalue.d_label = nodevalue_label;

   return;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void CommGraphWriter::writeGraphToTextStream(
   size_t record_number,
   std::ostream &os ) const
{
   /*
    * Gather graph data on d_root_rank and write out.
    */
   TBOX_ASSERT( record_number < d_records.size() );

   const Record &record = d_records[record_number];

   tbox::MessageStream ostr;

   for ( size_t inodev=0; inodev<record.d_node_values.size(); ++inodev ) {
      const NodeValue &nodev = record.d_node_values[inodev];
      ostr << nodev.d_value;
   }
   for ( size_t iedge=0; iedge<record.d_edges.size(); ++iedge ) {
      const Edge &edge = record.d_edges[iedge];
      ostr << edge.d_value << edge.d_dir << edge.d_other_node;
   }

   std::vector<char> tmpbuf( record.d_mpi.getRank() == d_root_rank ?
      ostr.getCurrentSize()*record.d_mpi.getSize() : 0 );

   if ( ostr.getCurrentSize() > 0 ) {
      record.d_mpi.Gather(
         (void*)ostr.getBufferStart(),
         int(ostr.getCurrentSize()),
         MPI_CHAR,
         (record.d_mpi.getRank() == d_root_rank ? &tmpbuf[0] : NULL),
         int(record.d_mpi.getRank() == d_root_rank ? ostr.getCurrentSize() : 0),
         MPI_CHAR,
         d_root_rank );
   }

   os.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
   os.precision(8);

   if ( record.d_mpi.getRank() == d_root_rank ) {

      os << "\nCommGraphWriter begin record number " << record_number << '\n';
      os << "# proc" << '\t' << "dir" << '\t' << "remote" << '\t' << "value" << '\t' << "label\n";

      if ( !tmpbuf.empty() ) {
         tbox::MessageStream istr( tmpbuf.size(),
                                   tbox::MessageStream::Read,
                                   &tmpbuf[0],
                                   false );

         for ( int src_rank = 0; src_rank < record.d_mpi.getSize(); ++src_rank ) {

            NodeValue tmpnodev;
            for ( size_t inodev=0; inodev<record.d_node_values.size(); ++inodev ) {
               istr >> tmpnodev.d_value;
               os << src_rank
                  << '\t' << tmpnodev.d_value
                  << '\t' << record.d_node_values[inodev].d_label
                  << '\n';
            }

            Edge tmpedge;
            for ( size_t iedge=0; iedge<record.d_edges.size(); ++iedge ) {
               istr >> tmpedge.d_value >> tmpedge.d_dir >> tmpedge.d_other_node;
               os << src_rank
                  << '\t' << (tmpedge.d_dir == FROM ? "<-" : "->")
                  << '\t' << tmpedge.d_other_node
                  << '\t' << tmpedge.d_value
                  << '\t' << record.d_edges[iedge].d_label
                  << '\n';
            }

         }
      }

      os << "CommGraphWriter end record number " << record_number << '\n';

   }

   return;
}



}
}
#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
