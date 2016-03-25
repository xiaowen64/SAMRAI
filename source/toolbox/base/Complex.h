//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/Complex.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	dcomplex class for old-style complex and new complex<double>
//

#ifndef included_tbox_Complex
#define included_tbox_Complex

#include "SAMRAI_config.h"

#ifndef included_complex
#include <complex>
#define included_complex
#endif

/**
 * @page toolbox_complex Toolbox Complex Type
 * 
 * @brief dcomplex is a typedef to overcome C++ compiler issues with
 * the std::complex type.
 *
 * The std::complex type should be a template however some older C++ compilers
 * implement complex as a double complex.  dcomplex is used to hide this
 * platform issue behind a typedef.
 *
 * NOTE: This should be removed when no longer required.
 *
 */

#ifndef LACKS_TEMPLATE_COMPLEX
typedef std::complex<double> dcomplex;
#else
typedef std::complex dcomplex;
#endif

#endif
