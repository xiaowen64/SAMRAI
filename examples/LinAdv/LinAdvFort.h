//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/LinAdv/LinAdvFort.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: F77 external declarations for SAMRAI linear advection example.
//

#include <math.h>
#include <signal.h>

extern "C" {
 
  void linadvinit_( 
  const int& , const double*, const double*,  const double*,
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& ,
#if (NDIM>1)
  const int& ,  
#endif
#if (NDIM>2)
  const int& ,  
#endif
  double*      , 
  const int&, 
  const double* , const double* );
  
  void linadvinitsine_(
  const int& , const double*, const double*,
  const double*, const double*,
  const int& , const int& ,
#if (NDIM>1)
  const int& , const int& ,
#endif
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int& ,
#if (NDIM>1)
  const int& ,
#endif
#if (NDIM>2)
  const int& ,
#endif
  double*      ,
  const int&,
  const double* , const double* ,
  const double&,  const double*);

  void initsphere_(
  const int& , const double*, const double*,  const double*,
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int&,
#if (NDIM>1)
  const int& ,  
#endif
#if (NDIM>2)
  const int& ,  
#endif
  double*      , 
  const double&, const double&, 
  const double*, const double&);

 
  void   stabledt_(
  const double*,
  const int& , const int& ,
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& ,
#if (NDIM>1)
  const int& ,
#endif
#if (NDIM>2)
  const int& ,
#endif
  const double*,
  const double*, 
  double&);
 
  void inittraceflux_(
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*, 
#if (NDIM>1)
  double*      , double*      , double*      , 
#endif
#if (NDIM>2)
  double*      , double*      , double*      , 
#endif
  double*      , double*      , double*      ); 
 
  void chartracing0_(
  const double&, const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& , const double&, const double&, const int& ,
  const double*,
  double*      , double*      ,
  double*      , double*      , 
  double*      , double*);
 
#if (NDIM>1)
  void chartracing1_(
  const double&, const int& , const int& , const int& ,const int& ,
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int& , const double&, const double&, const int& ,
  const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*);
 
#if (NDIM>2)
  void chartracing2_(
  const double&, const int& , const int& , const int& ,const int& ,
  const int& , const int& , 
  const int& , const double&, const double&, const int& ,
  const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*);
#endif
#endif
 
  void fluxcalculation_(
  const double&, const int& , const int& , 
#if (NDIM>2)
  const int& , 
#endif
  const double*,
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*,
  const double*, 
#if (NDIM>2)
  double*      , double*      , double*      , 
#endif
#if (NDIM>1)
  double*      , double*      , double*      , 
#endif
  double*      , double*      , double*      ); 
 
#if (NDIM>2)
  void fluxcorrec2d_(
  const double&, const int& , const int& , const int& ,const int& ,
  const int& , const int& , 
  const double*, const double*, const int&   ,
  const double*,
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  double*      , double*      , double*      , 
  double*      , double*      , double*      ); 

  void fluxcorrec3d_(
  const double&, const int& , const int& , const int& ,const int& ,
  const int& , const int& , 
  const double*, const double*, 
  const double*,
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  double*      , double*      , double*      , 
  double*      , double*      , double*      ); 
#endif

#if (NDIM==2)
  void fluxcorrec_(
  const double&, const int& , const int& , const int& ,const int& ,
  const double*,
  const double*, const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*      );
#endif

  void   consdiff_(
  const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*,
  const double*, const double*, 
#if (NDIM>1)
  const double*,
#endif
#if (NDIM>2)
  const double*,
#endif
  double*      );

  void getbdry_( const int& ,
  const int& , const int& , const int& , const int& , 
#if (NDIM>1)
  const int& , const int& , const int& , const int& , 
#endif
#if (NDIM>2)
  const int& , const int& , const int& , const int& , 
#endif
  const int& ,
#if (NDIM>1)
  const int& , 
#endif 
#if (NDIM>2)
  const int& , 
#endif 
  const int& ,
  const double*, const double&,
  double*      , 
  const double*, const double*, const int&);

#if (NDIM>2)
  void onethirdstate_(
  const double&, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, const double*, 
  const double*, const double*, const double*, 
  double*      ); 

  void fluxthird_(
  const double&, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, const double*, 
  const double*, 
  double*      , double*      , double*      );

  void fluxcorrecjt_(
  const double&, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, const double*,
  const double*, const double*, const double*,
  double*      , double*      , double*      ,
  double*      , double*      , double*      );
#endif

   void detectgrad_(
#if (NDIM == 2)
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#endif
#if (NDIM == 3)
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , 
      const int&, const int&,
      const double*,
      int* , int* );

   void detectshock_(
#if (NDIM == 2)
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#endif
#if (NDIM == 3)
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , const double& , 
      const int&, const int&,
      const double*,
      int* , int* );

   void stufprobc_(
      const int& , const int& , const int& , 
      const int& , const int& , const int& , const int& ,
      const int& , const int& , const int& );

}
