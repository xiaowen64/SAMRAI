/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/vtkcvis/vtkcvis.h $
  Language:  C++
  Date:      $LastChangedDate: 2005-12-20 13:26:41 -0800 (Tue, 20 Dec 2005) $
  Version:   $LastChangedRevision: 816 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkmyCommonWin32Header - manage Windows system differences
// .SECTION Description
// The vtkmyCommonWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkmyCommonWin32Header_h
#define __vtkmyCommonWin32Header_h

#include <vtkcvisConfigure.h>

#if defined(WIN32) && !defined(VTK_VTKCVIS_STATIC)
#if defined(vtkcvis_EXPORTS)
#define VTK_VTKCVIS_EXPORT __declspec( dllexport ) 
#else
#define VTK_VTKCVIS_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_VTKCVIS_EXPORT
#endif

#endif
