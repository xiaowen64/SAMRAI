/*
 * File:        PoissonSpecifications.C
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2005 The Regents of the University of California
 * Revision:    $Revision: 453 $
 * Modified:    $Date: 2005-06-16 10:19:28 -0700 (Thu, 16 Jun 2005) $
 * Description: Specifications for the scalar Poisson equation
 */



#include "PoissonSpecifications.h"

#ifdef DEBUG_NO_INLINE
#include "PoissonSpecifications.I"
#endif

namespace SAMRAI {
    namespace solv {




void PoissonSpecifications::printClassData( ostream &stream ) const
{
   stream << "PoissonSpecifications " << d_object_name << "\n"
          << "   D is ";
   if ( d_D_id != -1 ) {
      stream << "variable with patch id " << d_D_id << "\n";
   }
   else {
      stream << "constant with value " << d_D_constant << "\n";
   }
   stream << "   C is ";
   if ( d_C_zero ) {
      stream << "zero\n";
   }
   else if ( d_C_id != -1 ) {
      stream << "variable with patch id " << d_C_id << "\n";
   }
   else {
      stream << "constant with value " << d_C_constant << "\n";
   }
   return;
}


} // namespace SAMRAI
}
