/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dlbg/MDA_Access_instances.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: For explicit instantiations.
 */

#include "SAMRAI_config.h"
#include "MDA_Access.h"

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

template class MDA_Access< double, 2, MDA_OrderColMajor<2> >;
template class MDA_Access<    int, 2, MDA_OrderColMajor<2> >;

template class MDA_Access< double, 3, MDA_OrderColMajor<3> >;
template class MDA_Access<    int, 3, MDA_OrderColMajor<3> >;

#endif
