#ifndef included_tbox_ForAll
#define included_tbox_ForAll

namespace SAMRAI {
namespace tbox {

#if defined(HAVE_RAJA)

#include "RAJA/RAJA.hxx"

namespace parallel {
  typedef RAJA::cuda_exec<256> cuda;
}

template <typename policy, typename loop_body>
inline void for_all(int begin, int end, loop_body&& body) {
  RAJA::forall<policy>(begin, end, body);
}

#else

namespace parallel {
  typedef RAJA::seq_exec seq;
}

template <typename policy, typename loop_body>
inline void for_all(int begin, int end, loop_body&& body) {
  for(int i = begin; i < end; ++i) {
    body(i);
  }
}

#endif

}
}

#endif
