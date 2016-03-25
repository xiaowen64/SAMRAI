/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/solvers/poisson/PoissonSpecifications.C $
 * Package:     SAMRAI solvers
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Specifications for the scalar Poisson equation
 */



#include "PoissonSpecifications.h"

#ifdef DEBUG_NO_INLINE
#include "PoissonSpecifications.I"
#endif

namespace SAMRAI {
    namespace solv {




void PoissonSpecifications::printClassData( std::ostream &stream ) const
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
