//
// File:	$URL$
// Package:	SAMRAI templates
// Copyright:	(c) 1997-2004 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2207 $
// Modified:	$LastChangedDate: 2008-06-05 15:21:40 -0700 (Thu, 05 Jun 2008) $
// Description:	Empty class that is built so dimensional archives have
//              something in them.  If implicit templating is used then
//              some of the packages have no objects.   This causes
//              problems on some machines.
//


namespace SAMRAI {
   namespace mesh {
      template<int DIM> class Empty {
      };

      template class Empty< NDIM >;
   }
}


