/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple utility class for interfacing with MPI 
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const SAMRAI_MPI& SAMRAI_MPI::getSAMRAIWorld() {
   return s_samrai_world;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
SAMRAI_MPI::SAMRAI_MPI(
   const Comm& comm):
   d_comm(comm),
   d_rank(0),
   d_size(1)
{
#ifdef HAVE_MPI
   if (comm != MPI_COMM_NULL) {
      if (s_mpi_is_initialized) {
         MPI_Comm_rank(d_comm, &d_rank);
         MPI_Comm_size(d_comm, &d_size);
      }
   }
#endif
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void SAMRAI_MPI::dupCommunicator(
   const SAMRAI_MPI& r)
{
#ifdef HAVE_MPI
   int rval = r.Comm_dup(&d_comm);
   if (rval != MPI_SUCCESS) {
      TBOX_ERROR("SAMRAI_MPI::dupCommunicator: Error duplicating\n"
         << "communicator.");
   }
   MPI_Comm_rank(d_comm, &d_rank);
   MPI_Comm_size(d_comm, &d_size);
   TBOX_ASSERT(d_rank == r.d_rank);
   TBOX_ASSERT(d_size == r.d_size);
#else
   NULL_USE(r);
#endif
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void SAMRAI_MPI::freeCommunicator()
{
#ifdef HAVE_MPI
   if ( d_comm != MPI_COMM_NULL ) {
      TBOX_ASSERT( tbox::SAMRAI_MPI::usingMPI() );
      Comm_free(&d_comm);
      // d_comm is now set to MPI_COMM_NULL;
   }
#else
   d_comm = MPI_COMM_NULL;
#endif
   d_rank = 0;
   d_size = 1;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const SAMRAI_MPI& SAMRAI_MPI::operator = (
   const SAMRAI_MPI& rhs)
{
   d_comm = rhs.d_comm;
   d_rank = rhs.d_rank;
   d_size = rhs.d_size;
   return *this;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool SAMRAI_MPI::operator == (
   const SAMRAI_MPI& rhs) const
{
   return d_comm == rhs.d_comm;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool SAMRAI_MPI::operator != (
   const SAMRAI_MPI& rhs) const
{
   return d_comm != rhs.d_comm;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int SAMRAI_MPI::getRank() const
{
   return d_rank;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int SAMRAI_MPI::getSize() const
{
   return d_size;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
int SAMRAI_MPI::getInvalidRank()
{
   return s_invalid_rank;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const SAMRAI_MPI::Comm& SAMRAI_MPI::getCommunicator() const
{
   return d_comm;
}

/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void SAMRAI_MPI::setCommunicator(
   const Comm& comm)
{
   d_comm = comm;
   d_rank = 0;
   d_size = 1;
#ifdef HAVE_MPI
   if (s_mpi_is_initialized) {
      if (d_comm != MPI_COMM_NULL) {
         MPI_Comm_rank(d_comm, &d_rank);
         MPI_Comm_size(d_comm, &d_size);
      }
   }
#endif
}

#ifndef HAVE_MPI
/*
 **************************************************************************
 **************************************************************************
 */
SAMRAI_INLINE_KEYWORD
SAMRAI_MPI::Status::Status():
   MPI_SOURCE(-1),
   MPI_TAG(-1),
   MPI_ERROR(-1)
{
}
#endif

}
}
