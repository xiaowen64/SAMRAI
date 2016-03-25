//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/base/Complex.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1735 $
// Modified:	$LastChangedDate: 2007-12-05 15:01:59 -0800 (Wed, 05 Dec 2007) $
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
