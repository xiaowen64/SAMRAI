#include "SAMRAI_config.h"
#include "MDA_Access.h"

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

template class MDA_Access< double, 2, MDA_OrderColMajor<2> >;
template class MDA_Access<    int, 2, MDA_OrderColMajor<2> >;

template class MDA_Access< double, 3, MDA_OrderColMajor<3> >;
template class MDA_Access<    int, 3, MDA_OrderColMajor<3> >;

#endif
