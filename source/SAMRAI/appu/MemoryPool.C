#include "MemoryPool.h"

#if defined(HAVE_CUDA)

namespace SAMRAI {
namespace appu {

MemoryPool * MemoryPool::s_memory_pool_instance(0);

MemoryPool *
MemoryPool::getMemoryPool()
{
  if (!s_memory_pool_instance) {
    s_memory_pool_instance = new MemoryPool();
  }
  return s_memory_pool_instance;
}

MemoryPool::MemoryPool()
{
  d_cnmem_device.device = 0;
  d_cnmem_device.size = 32 * 1024 * 1024;
  d_cnmem_device.numStreams = 0;
  d_cnmem_device.streams = NULL;
  d_cnmem_device.streamSizes = NULL;

  cnmemInit(1, &d_cnmem_device, CNMEM_FLAGS_DEFAULT);
}

MemoryPool::~MemoryPool()
{
  cnmemFinalize();
}

void * MemoryPool::alloc(size_t size)
{
  void * ret;
  cnmemMalloc(&ret, size, 0);
  return ret;
}

void MemoryPool::free(void * ptr)
{
  cnmemFree(ptr, 0);
}


}
}

#endif
