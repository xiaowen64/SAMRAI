/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   STL allocator for allocating unified memory when available
 *
 ************************************************************************/

#ifndef included_pdat_ManagedAllocator
#define included_pdat_ManagedAllocator

#include "SAMRAI/SAMRAI_config.h"

#include <cstdlib>
#include <cstddef>
#if defined(HAVE_CUDA)
#include <cuda_runtime_api.h>
#endif

#if defined(ENABLE_SIMPOOL)
#include "DynamicPoolAllocator.hpp"
#endif

namespace SAMRAI {
namespace pdat {

struct UMAllocator
{
#if defined(HAVE_CUDA)
  static inline void *allocate(std::size_t size) { void *ptr; cudaMallocManaged(&ptr, size); return ptr; }
  static inline void deallocate(void *ptr) { cudaFree(ptr); }
#else
  static inline void *allocate(std::size_t size) { return std::malloc(size); }
  static inline void deallocate(void *ptr) { std::free(ptr); }
#endif
};

template <class T>
struct ManagedAllocator {
  typedef T value_type;
  typedef std::size_t size_type;

#if defined(ENABLE_SIMPOOL)
  typedef DynamicPoolAllocator<UMAllocator> PoolType;
  PoolType &m;
  ManagedAllocator() : m(PoolType::getInstance()) { }
  ManagedAllocator(const ManagedAllocator& a) : m(a.m) { }
#else
  ManagedAllocator() {}
  ManagedAllocator(const ManagedAllocator& a) { }
#endif

  T* allocate(std::size_t n) {
#if defined(ENABLE_SIMPOOL)
    return static_cast<T*>( m.allocate( n * sizeof(T) ) );
#else
    return static_cast<T*>( UMAllocator::allocate( n * sizeof(T) ) );
#endif
  }

  void deallocate(T* p, std::size_t n) {
#if defined(ENABLE_SIMPOOL)
    m.deallocate(p);
#else
    UMAllocator::deallocate(p);
#endif
  }

  size_type max_size() const { return std::numeric_limits<size_type>::max(); }
};

template <class T, class U>
bool operator==(const ManagedAllocator<T>&, const ManagedAllocator<U>&) { return true; }

template <class T, class U>
bool operator!=(const ManagedAllocator<T>&, const ManagedAllocator<U>&) { return false; }

}
}

#endif
