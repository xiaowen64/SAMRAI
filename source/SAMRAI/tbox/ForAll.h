#ifndef included_tbox_ForAll
#define included_tbox_ForAll

#if defined(HAVE_RAJA)

#include "RAJA/RAJA.hpp"

// namespace SAMRAI {
// namespace tbox {
// 
// namespace parallel {
//   typedef RAJA::cuda_exec<256> cuda;
// }
// 
// template <typename policy, typename loop_body>
// inline void for_all(int begin, int end, loop_body&& body) {
//   RAJA::forall<policy>(begin, end, body);
// }
// 
// #else
// 
// namespace SAMRAI {
// namespace tbox {
// 
// namespace parallel {
//   typedef RAJA::seq_exec seq;
// }
// 
// template <typename policy, typename loop_body>
// inline void for_all(int begin, int end, loop_body&& body) {
//   for(int i = begin; i < end; ++i) {
//     body(i);
//   }
// }
// 
// #endif
// 
// }
// }

#endif
#endif
