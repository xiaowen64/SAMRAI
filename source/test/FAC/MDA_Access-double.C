/*
 * File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/MDA_Access-double.C $
 * Package:     SAMRAI tests
 * Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
 * Revision:    $LastChangedRevision: 1704 $
 * Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
 * Description: Explicit instantiations in FAC solver test.
 */

#include "SAMRAI_config.h"
#include "MDA_Access.h"

template class MDA_Access     < double, 2, MDA_OrderColMajor<2> >;
template class MDA_AccessConst< double, 2, MDA_OrderColMajor<2> >;

