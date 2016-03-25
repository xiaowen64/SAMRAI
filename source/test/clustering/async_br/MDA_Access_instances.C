/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/clustering/async_br/MDA_Access_instances.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
