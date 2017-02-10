#ifndef included_appu_MemoryPool
#define included_appu_MemoryPool

#include "SAMRAI/SAMRAI_config.h"

#if defined(HAVE_CUDA) && defined(HAVE_CNMEM)


#include "cnmem.h"

namespace SAMRAI {
namespace appu {

class MemoryPool
{
  public:
    static MemoryPool * getMemoryPool();

    void* alloc(size_t size);

    void free(void * ptr);

  protected:
    MemoryPool();

    virtual ~MemoryPool();

  private:
    static MemoryPool* s_memory_pool_instance;
    cnmemDevice_t d_cnmem_device;
};

}
}

#endif
#endif
