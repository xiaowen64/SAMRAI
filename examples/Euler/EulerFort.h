//
// File:        EulerFort.h
// Package:     SAMRAI application
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
