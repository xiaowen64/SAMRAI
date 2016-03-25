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
// .NAME vtkBar - Bar class for vtk
// .SECTION Description
// None.

#ifndef __vtkBar_h
#define __vtkBar_h

#include <vtkObject.h>
#include "vtkcvis.h"

class VTK_VTKCVIS_EXPORT vtkBar : public vtkObject
{
public:
  static vtkBar *New();
  vtkTypeRevisionMacro(vtkBar,vtkObject);

protected:
  vtkBar() {};
  ~vtkBar() {};
private:
  vtkBar(const vtkBar&);  // Not implemented.
  void operator=(const vtkBar&);  // Not implemented.
};

#endif 
