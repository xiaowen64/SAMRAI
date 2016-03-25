/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/FAC/MDA_Access-int.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1917 $
 * Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
 * Description: Explicit instantiations in FAC solver test.
 */

#include "SAMRAI_config.h"
#include "MDA_Access.h"

template class MDA_Access     < int, 2, MDA_OrderColMajor<2> >;
template class MDA_AccessConst< int, 2, MDA_OrderColMajor<2> >;

