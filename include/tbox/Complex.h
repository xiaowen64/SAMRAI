//
// File:	Complex.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	dcomplex class for old-style complex and new complex<double>
//

#ifndef included_tbox_Complex
#define included_tbox_Complex

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_complex
#include <complex>
#define included_complex
using namespace std;
#endif

#ifndef LACKS_TEMPLATE_COMPLEX
typedef complex<double> dcomplex;
#else
typedef complex dcomplex;
#endif

#endif
