/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/vtkcvis/vtkBar.h $
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
