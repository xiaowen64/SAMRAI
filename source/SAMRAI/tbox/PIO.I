/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Parallel I/O classes pout, perr, and plog and control class
 *
 ************************************************************************/

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * Close the output streams.  Flush both cout and cerr.  If logging,
 * then flush and close the log stream.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void PIO::finalize()
{
   std::cout.flush();
   std::cerr.flush();
   shutdownFilestream();
}

}
}
