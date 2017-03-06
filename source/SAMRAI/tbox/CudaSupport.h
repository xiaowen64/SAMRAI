#ifndef included_tbox_CudaSupport
#define included_tbox_CudaSupport

#define SAMRAI_INLINE inline

#if defined(HAVE_CUDA)

#if defined(__CUDACC__)

#define SAMRAI_DEVICE __device__
#define SAMRAI_HOST_DEVICE __host__ __device__

#else

#define SAMRAI_DEVICE 
#define SAMRAI_HOST_DEVICE 

#endif

#endif
#endif
