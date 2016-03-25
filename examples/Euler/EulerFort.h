//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/Euler/EulerFort.h $
// Package:     SAMRAI application
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: F77 external declarations for SAMRAI Euler gas dynamics ex.
//

extern "C" {

  void eulerinit_(
  const int& , const double*, const double*, const double*,
  const int& , const int& ,
  const int& , const int& ,
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int&,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  const double&,
  double*      , double*      , double*      ,
  const int&,
  const double*,
  const double*, const double*, const double*);

  void eulerinitsphere_(
  const int& , const double*, const double*, const double*,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int&,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  const double&,
  double*      , double*      , double*      ,
  const double&, const double*, const double&,
  const double&, const double*, const double&,
  const double*, const double&);

  void   stabledt_(
  const double*,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const int& ,
  const int& ,
#if (NDIM>2)
  const int& ,
#endif
  const double&,
  const double*, const double*, const double*, double&);
 
  void inittraceflux_(
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*, const double*, const double*,
  double*      , double*      , double*      , 
#if (NDIM>2)
  double*      , double*      , double*      , 
#endif
  double*      , double*      , double*      );
 
  void computesound_(
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double&, 
  const double*, const double*, const double*,
  double*                      ); 
 
  void chartracing0_(
  const double&,
  const int& , const int& ,
  const int& , const int& ,
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int&    , const double& , const double& , const int&    ,
  const double* ,
  double*       , double*       ,
  double*       , double*       ,
  double*       , 
  double*       , double*       );

  void chartracing1_(
  const double&, 
  const int& , const int& , 
  const int& ,const int& ,
#if (NDIM>2)
  const int& , const int& ,
#endif
  const int&   , const double&, const double&, const int&    ,
  const double*,
  double*      , double*      ,
  double*      , double*      , 
  double*      ,  
  double*      , double*);

#if (NDIM == 3)
  void chartracing2_(
  const double&, 
  const int& , const int& , 
  const int& , const int& ,
  const int& , const int& ,
  const int& ,   const double&, const double&, const int&    ,
  const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      ,
  double*      , double*);
#endif
 
  void fluxcalculation_(
  const double&, const int& ,
#if (NDIM>2)
  const int& , 
#endif
  const int& ,
  const double*,
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double&,
  const int&,
  const double*, const double*, const double*,
#if (NDIM>2)
  double*      , double*      , double*      , 
#endif
  double*      , double*      , double*      , 
  double*      , double*      , double*      ); 
 
#if (NDIM == 3)
  void fluxcorrec2d_(
  const double&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, const double&, const int&,
  const double*,
  const double*, 
  const double*,
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  double*      , double*      , double*      , 
  double*      , double*      , double*      );

  void fluxcorrec3d_(
  const double&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double*, const double&, 
  const double*, const double*, const double*,
  const double*, const double*, const double*, 
  const double*, const double*, const double*, 
  double*      , double*      , double*      , 
  double*      , double*      , double*      ); 
#endif

#if (NDIM==2)
  void fluxcorrec_(
  const double&, 
  const int& , const int& , const int& ,const int& ,
  const double*, const double&,
  const double*, const double*, const double*,
  double*      , double*      ,
  double*      , double*      ,
  double*      , double*      );
#endif

  void   consdiff_(
  const int& , const int& , 
  const int& , const int& , 
#if (NDIM>2)
  const int& , const int& , 
#endif
  const double*,
  const double*,
  const double*,
#if (NDIM>2)
  const double*,
#endif
  const double&, 
  double*      , double*      , double*      );

#if (NDIM>2)
  void onethirdstate_(
  const double&, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double&, 
  const double*, const double*, const double*,
  const double*, const double*, const double*,
  double*      ); 

  void fluxthird_(
  const double&, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double&, 
  const int&,
  const double*, const double*, const double*, const double*, 
  double*      , double*,       double*);

  void fluxcorrecjt_(
  const double&, const double*, const int&,
  const int& , const int& , const int& ,const int& , const int& , const int& , 
  const double&, 
  const double*, const double*, const double*,
  const double*, const double*, const double*,
  double*      , double*      , double*      ,
  double*      , double*      , double*      );
#endif

#if (NDIM == 2)
   void conservlinint2d_(
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*,
      double*, double*, double*, const int&,
      double*, double*, double*,
      double*, double*,
      double*, double*, double*, double*);

   void conservavg2d_(
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*);
#endif
#if (NDIM == 3)
   void conservlinint3d_(
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*,
      double*, double*, double*, const int&,
      double*, double*, double*,
      double*, double*, double*,
      double*, double*, double*, double*, double*, double*);

   void conservavg3d_(
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int*, const double*, const double*, const double&,
      const double*, const double*,
      const double*, const double*,
      double*, double*,
      double*);
#endif

   void detectgrad_(
      const int& , const int& , 
      const int& , const int& , 
#if (NDIM>2)
      const int& , const int& , 
#endif
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#if (NDIM>2)
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , 
      const int&, const int&,
      const double*,
      int* , int* );

   void detectshock_(
      const int& , const int& , 
      const int& , const int& , 
#if (NDIM>2)
      const int& , const int& , 
#endif
      const int& , const int& , const int&, 
      const int& , const int& , const int&, 
#if (NDIM>2)
      const int& , const int& , const int&, 
#endif
      const double* , 
      const double& , const double& , 
      const int&, const int&,
      const double*,
      int* , int* );

   void stufprobc_(
      const int& , const int& , const int& , 
      const int& , const int& , const int& , 
      const int& , const int& ,
      const int& , const int& , const int& );
}
