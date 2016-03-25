/*=========================================================================

  Program:   Vizamrai
  Date:      $Date: 2004-11-30 08:48:55 -0800 (Tue, 30 Nov 2004) $
  Version:   $Revision: 5 $

  Based on the hedgehog code from VTK with some local modifications.

  Program:   Visualization Toolkit
  Module:    RCSfile: vtkHedgeHog.h,v 
  Language:  C++
  Date:      Date: 2002/01/22 15:29:23 
  Version:   Revision: 1.34 

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkcvisHedgeHog - create oriented lines from vector data
// .SECTION Description
// vtkcvisHedgeHog creates oriented lines from the input data set. Line
// length is controlled by vector (or normal) magnitude times scale
// factor. Lines are colored by scalar data, if available.

#ifndef __vtkcvisHedgeHog_h
#define __vtkcvisHedgeHog_h

#include "vtkDataSetToPolyDataFilter.h"

#define VTK_USE_VECTOR 0
#define VTK_USE_NORMAL 1

class VTK_GRAPHICS_EXPORT vtkcvisHedgeHog : public vtkDataSetToPolyDataFilter
{
public:
  static vtkcvisHedgeHog *New();
  vtkTypeRevisionMacro(vtkcvisHedgeHog,vtkDataSetToPolyDataFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set scale factor to control size of oriented lines.
  vtkSetVector2Macro(VectorScaleFactor,float);
  vtkGetVector2Macro(VectorScaleFactor,float);

  // Description:
  // Set the range to scale to, needed since we want same
  // scaling to occur over multiple boxes.
  vtkSetVector2Macro(ScalingRange,float);
  vtkGetVector2Macro(ScalingRange,float);

  // Description:
  // Control cropping
  vtkSetMacro(Crop,int);
  vtkGetMacro(Crop,int);

  // Description:
  // Set cropping min/max values for velocities
  vtkSetVector2Macro(CroppingRange,float);
  vtkGetVector2Macro(CroppingRange,float);

protected:
  vtkcvisHedgeHog();
  ~vtkcvisHedgeHog() {};

  void Execute();
  float VectorScaleFactor[2];
  float ScalingRange[2];
  int Crop;
  float CroppingRange[2];

private:
  vtkcvisHedgeHog(const vtkcvisHedgeHog&);  // Not implemented.
  void operator=(const vtkcvisHedgeHog&);  // Not implemented.
};

#endif



