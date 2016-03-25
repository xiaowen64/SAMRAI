/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2004-11-30 08:48:55 -0800 (Tue, 30 Nov 2004) $
  Version:   $Revision: 5 $

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
