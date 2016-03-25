//
// File:        $RCSfile$
// Package:     SAMRAI application
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Fortran function declarations for modified-Bratu problem
//

// Header file to define interfaces between Fortran and C++ code for
// the ModifiedBratuProblem class. 

#ifndef included_modifiedBratuFort
#define included_modifiedBratuFort

#define FORT_FILL            fill_
#define FORT_EVALBRATU       evalbratu_
#define FORT_EVALDIFFUSION   evaldiffusioncoef_
#define FORT_EVALEXPONENTIAL evalexponential_
#define FORT_EVALFACEFLUXES  evalfacefluxes_
#define FORT_EVALSOURCE      evalsource_
#define FORT_EWBCFLUXFIX     ewbcfluxfix_
#define FORT_NSBCFLUXFIX     nsbcfluxfix_
#define FORT_TBBCFLUXFIX     tbbcfluxfix_
#define FORT_EWFLUXCOPY      ewfluxcopy_
#define FORT_NSFLUXCOPY      nsfluxcopy_
#define FORT_TBFLUXCOPY      tbfluxcopy_
#define FORT_BRATUJV         bratujv_
#define FORT_SETBC           setbc_
#define FORT_ERROR           error_
#define FORT_EVALF           evalf_
#define FORT_PROLONG         prolong_

extern "C"
{

#if (NDIM == 1)

   void FORT_FILL( const int&, const int&, 
		   const double*,
		   const double*, const double*, const double*,
		   const int& );

   void FORT_EVALBRATU( const int&, const int&, 
			const int&,
			const double*, 
			const double*, const double*,
			const double*,
			const double*,
			const double*, const double&,
			const double* );

   void FORT_EVALDIFFUSION( const int&, const int&, 
			    const double*, const double*, const double*,
			    const double* );

   void FORT_EVALEXPONENTIAL( const int&, const int&, 
			      const double*,
			      const double&,
			      const double* );

   void FORT_EVALFACEFLUXES( const int&, const int&, 
			     const int&,
			     const double*,
			     const double*,
			     const double*,
			     const double* );

   void FORT_EVALSOURCE( const int&, const int&, 
			 const double&,
			 const double*, const double*, const double*,
			 const double&,
			 const double* );

   void FORT_EWBCFLUXFIX( const int&, const int&, 
			  const int&,
			  const double*,
			  const double*,
			  const double*,
			  const int*, const int*,
			  const int& );

   void FORT_EWFLUXCOPY( const int&, const int&, 
			 const double*,
			 const double*,
			 const int& );

   void FORT_BRATUJV( const int&, const int&, 
		      const int&,
		      const double*,
		      const double*,
		      const double*,
		      const double*,
		      const double&,
		      const double* );

   void FORT_SETBC( const int&, const int&, 
		    const int&,
		    const double*,
		    const int*, const int*,
		    const int& );

   void FORT_ERROR( const int&, const int&, 
		    const double*, const double*,
		    const double&,
		    const double*, const double*, const double*,
		    const double&,
		    const double&, 
		    const double& );

   void FORT_EVALF( const int&, const int&, 
		    const double*, 
		    const double*,
		    const double*,
		    const double* );

   void FORT_VAXPY( const int&, const int&, 
		    const double*, 
		    const double*,
		    const double* );

   void compfacdiag_(const int& , const int& ,
		     const double& , const double& ,
		     const double* , const double* ,
		     double*);

   void compfacoffdiag_(const int& , const int& ,
			const double& , const double& ,
			const double* , 
			double* );

   void compdiffcoef_(const int& , const int& ,
		      const double* , const double* ,
		      const double& ,
		      double*);

   void compexpu_(const int& , const int& ,
		  const int& ,
		  const double& ,
		  const double* ,
		  double*);

   void compsrc_(const int& , const int& ,
		 const int& ,
		 const double* , const double* , const double& ,
		 const double& , 
		 const double* ,
		 double*);

   void compsrcderv_(const int& , const int& ,
		     const int& ,
		     const double* , const double* , const double& ,
		     const double& , 
		     const double* ,
		     double*);

   void compsideflux_(const int& , const int& ,
		      const int& ,
		      const double* , 
		      const double* ,
		      const double* ,
		      double*);

   void fluxbdryfix_(const int& , const int& ,
		     const int& , const int& ,
		     const int& ,
		     const int& ,
		     const double& ,
		     double*);

   void fluxcopy0_(const int& , const int& ,
		   const int& ,
		   const double* ,
		   double*);

   void compresidual_(const int& , const int& ,
		      const int& ,
		      const double* , const double& ,
		      const double* ,
		      const double* ,
		      const double* ,
		      const double* ,
		      const double* , 
		      double*);


// Bonus function

   void adjcrsfineoffdiag1d_( const int&, const int&, 
			      const int&, const int&, 
			      const int&, const int&,
			      double* );
#endif
#if (NDIM == 2)

   void FORT_FILL( const int&, const int&, 
		   const int&, const int&, 
		   const double*,
		   const double*, const double*, const double*,
		   const int& );

   void FORT_EVALBRATU( const int&, const int&, 
			const int&, const int&, 
			const int&,
			const double*, const double*,
			const double*, const double*,
			const double*,
			const double*,
			const double*, const double&,
			const double* );

   void FORT_EVALDIFFUSION( const int&, const int&, 
			    const int&, const int&, 
			    const double*, const double*, const double*,
			    const double*, const double* );

   void FORT_EVALEXPONENTIAL( const int&, const int&, 
			      const int&, const int&, 
			      const double*,
			      const double&,
			      const double* );

   void FORT_EVALFACEFLUXES( const int&, const int&, 
			     const int&, const int&, 
			     const int&,
			     const double*, const double*,
			     const double*,
			     const double*,
			     const double*, const double* );

   void FORT_EVALSOURCE( const int&, const int&, 
			 const int&, const int&, 
			 const double&,
			 const double*, const double*, const double*,
			 const double&,
			 const double* );

   void FORT_EWBCFLUXFIX( const int&, const int&, 
			  const int&, const int&, 
			  const int&,
			  const double*,
			  const double*,
			  const double*,
			  const int*, const int*,
			  const int& );

   void FORT_NSBCFLUXFIX( const int&, const int&, 
			  const int&, const int&,
			  const int&,
			  const double*,
			  const double*,
			  const double*,
			  const int*, const int*,
			  const int& );

   void FORT_EWFLUXCOPY( const int&, const int&, 
			 const int&, const int&,
			 const double*,
			 const double*,
			 const int& );
	       
   void FORT_NSFLUXCOPY( const int&, const int&, 
			 const int&, const int&, 
			 const double*,
			 const double*,
			 const int& );

   void FORT_BRATUJV( const int&, const int&, 
		      const int&, const int&,
		      const int&,
		      const double*, const double*,
		      const double*,
		      const double*,
		      const double*,
		      const double&,
		      const double* );

   void FORT_SETBC( const int&, const int&, 
		    const int&, const int&,
		    const int&,
		    const double*,
		    const int*, const int*,
		    const int& );

   void FORT_ERROR( const int&, const int&, 
		    const int&, const int&,
		    const double*, const double*,
		    const double&,
		    const double*, const double*, const double*,
		    const double&,
		    const double&, 
		    const double& );

   void FORT_EVALF( const int&, const int&, 
		    const int&, const int&,
		    const double*, 
		    const double*,
		    const double*,
		    const double* );

   void FORT_VAXPY( const int&, const int&, 
		    const int&, const int&,
		    const double*, 
		    const double*,
		    const double* );

   /* These functions are in FACjacobian.m4 */

   void compjv_( const int &ifirst0, const int &ilast0,
		 const int &ifirst1, const int &ilast1,
		 const int &gwc,
		 const double *diag,
		 const double *flux0, const double *flux1,
		 const double *v,
		 const double *dx,
		 const double &dt,
		 double *jv );

   void compfacdiag_(const int& , const int& ,
		     const int& , const int& ,
		     const double& , const double& ,
		     const double* , const double* ,
		     double*);

   void compfacoffdiag_(const int& , const int& ,
			const int& , const int& ,
			const double& , const double& ,
			const double* , const double* , 
			double* , double* );

   void compdiffcoef_(const int& , const int& ,
		      const int& , const int& ,
		      const double* , const double* ,
		      const double& ,
		      double* , double*);

   void compexpu_(const int& , const int& ,
		  const int& , const int& ,
		  const int& ,
		  const double& ,
		  const double* ,
		  double*);

   void compsrc_(const int& , const int& ,
		 const int& , const int& ,
		 const int& ,
		 const double* , const double* , const double& ,
		 const double& , 
		 const double* ,
		 double*);

   void compsrcderv_(const int& , const int& ,
		     const int& , const int& ,
		     const int& ,
		     const double* , const double* , const double& ,
		     const double& , 
		     const double* ,
		     double*);

   void compsideflux_(const int& , const int& ,
		      const int& , const int& ,
		      const int& ,
		      const double* , 
		      const double* , const double* ,
		      const double* ,
		      double* , double*);

   void fluxbdryfix_(const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& ,
		     const int& ,
		     const double& ,
		     double* , double*);

   void fluxcopy0_(const int& , const int& ,
		   const int& , const int& ,
		   const int& ,
		   const double* ,
		   double*);

   void fluxcopy1_(const int& , const int& ,
		   const int& , const int& ,
		   const int& ,
		   const double* ,
		   double*);

   void compresidual_(const int& , const int& ,
		      const int& , const int& ,
		      const int& ,
		      const double* , const double& ,
		      const double* ,
		      const double* ,
		      const double* ,
		      const double* ,
		      const double* , const double* , 
		      double*);
#endif
#if (NDIM == 3)
   void FORT_FILL( const int&, const int&, 
		   const int&, const int&, 
		   const int&,
		   const double*,
		   const double*, const double*, const double*,
		   const int& );

   void FORT_EVALBRATU( const int&, const int&, 
			const int&, const int&, 
			const int&, const int&,
			const int&,
			const double*, const double*, const double*,
			const double*, const double*,
			const double*,
			const double*,
			const double*, const double&,
			const double* );

   void FORT_EVALDIFFUSION( const int&, const int&, 
			    const int&, const int&, 
			    const int&, const int&,
			    const double*, const double*, const double*,
			    const double*, const double*, const double* );

   void FORT_EVALEXPONENTIAL( const int&, const int&, 
			      const int&, const int&, 
			      const int&, const int&,
			      const double*,
			      const double&,
			      const double* );

   void FORT_EVALFACEFLUXES( const int&, const int&, 
			     const int&, const int&, 
			     const int&, const int&,
			     const int&,
			     const double*, const double*, const double*,
			     const double*,
			     const double*,
			     const double*, const double*, const double* );

   void FORT_EVALSOURCE( const int&, const int&, 
			 const int&, const int&, 
			 const int&, const int&,
			 const double&,
			 const double*, const double*, const double*,
			 const double&,
			 const double* );

   void FORT_EWBCFLUXFIX( const int&, const int&, 
			  const int&, const int&, 
			  const int&, const int&,
			  const int&, 
			  const double*,
			  const double*,
			  const double*,
			  const int*, const int*,
			  const int& );

   void FORT_NSBCFLUXFIX( const int&, const int&, 
			  const int&, const int&, 
			  const int&, const int&,
			  const int&, 
			  const double*,
			  const double*,
			  const double*,
			  const int*, const int*,
			  const int& );

   void FORT_TBBCFLUXFIX( const int&, const int&, 
			  const int&, const int&, 
			  const int&, const int&,
			  const int&, 
			  const double*,
			  const double*,
			  const double*,
			  const int*, const int*,
			  const int& );

   void FORT_EWFLUXCOPY( const int&, const int&, 
			 const int&, const int&, 
			 const int&, const int&,
			 const double*,
			 const double*,
			 const int& );
	       
   void FORT_NSFLUXCOPY( const int&, const int&, 
			 const int&, const int&, 
			 const int&, const int&,
			 const double*,
			 const double*,
			 const int& );
	       
   void FORT_TBFLUXCOPY( const int&, const int&, 
			 const int&, const int&, 
			 const int&, const int&,
			 const double*,
			 const double*,
			 const int& );

   void FORT_BRATUJV( const int&, const int&, 
		      const int&, const int&, 
		      const int&, const int&,
		      const int&,
		      const double*, const double*, const double*,
		      const double*,
		      const double*,
		      const double*,
		      const double&,
		      const double* );

   void FORT_SETBC( const int&, const int&, 
		    const int&, const int&, 
		    const int&, const int&,
		    const int&,
		    const double*,
		    const int*, const int*,
		    const int& );

   void FORT_ERROR( const int&, const int&, 
		    const int&, const int&, 
		    const int&, const int&,
		    const double*, const double*,
		    const double&,
		    const double*, const double*, const double*,
		    const double&,
		    const double&, 
		    const double& );                

   void FORT_EVALF( const int&, const int&, 
		    const int&, const int&, 
		    const int&, const int&, 
		    const double*, 
		    const double*,
		    const double*,
		    const double* );

   void FORT_VAXPY( const int&, const int&, 
		    const int&, const int&, 
		    const int&, const int&, 
		    const double*, 
		    const double*,
		    const double* );

   /* These functions are in FACjacobian.m4 */

   void compjv_( const int &ifirst0, const int &ilast0,
		 const int &ifirst1, const int &ilast1,
		 const int &ifirst2, const int &ilast2,
		 const int &gwc,
		 const double *diag,
		 const double *flux0, const double *flux1, const double *flux2,
		 const double *v,
		 const double *dx,
		 const double &dt,
		 double *jv );

   void compfacdiag_(const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const double& , const double& ,
		     const double* , const double* ,
		     double*);
   void compfacoffdiag_(const int& , const int& ,
			const int& , const int& ,
			const int& , const int& ,
			const double& , const double& ,
			const double* , const double* , const double* ,
			double* , double* , double* );

   void compdiffcoef_(const int& , const int& ,
		      const int& , const int& ,
		      const int& , const int& ,
		      const double* , const double* ,
		      const double& ,
		      double* , double* , double*);

   void compexpu_(const int& , const int& ,
		  const int& , const int& ,
		  const int& , const int& ,
		  const int& ,
		  const double& ,
		  const double* ,
		  double*);               

   void compsrc_(const int& , const int& ,
		 const int& , const int& ,
		 const int& , const int& ,
		 const int& ,
		 const double* , const double* , const double& ,
		 const double& , 
		 const double* ,
		 double*);

   void compsrcderv_(const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& ,
		     const double* , const double* , const double& ,
		     const double& , 
		     const double* ,
		     double*);

   void compsideflux_(const int& , const int& ,
		      const int& , const int& ,
		      const int& , const int& ,
		      const int& ,
		      const double* , 
		      const double* , const double* , const double* ,
		      const double* ,
		      double* , double* , double*);

   void fluxbdryfix_(const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& , const int& ,
		     const int& ,
		     const int& ,
		     const double& ,
		     double* , double* , double*);

   void fluxcopy0_(const int& , const int& ,
		   const int& , const int& ,
		   const int& , const int& ,
		   const int& ,
		   const double* ,
		   double*);
   void fluxcopy1_(const int& , const int& ,
		   const int& , const int& ,
		   const int& , const int& ,
		   const int& ,
		   const double* ,
		   double*);

   void compresidual_(const int& , const int& ,
		      const int& , const int& ,
		      const int& , const int& ,
		      const int& ,
		      const double* , const double& ,
		      const double* ,
		      const double* ,
		      const double* ,
		      const double* ,
		      const double* , const double* , const double* ,
		      double*);
#endif

};

#endif
