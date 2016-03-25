/*=========================================================================
  Program:   Vizamrai
  Date:      $LastChangedDate: 2005-12-20 13:26:41 -0800 (Tue, 20 Dec 2005) $
  Version:   $LastChangedRevision: 816 $

  Based on the hedgehog code from VTK with some local modifications.

  Program:   Visualization Toolkit
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
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkcvisHedgeHog.h"
#include "vtkObjectFactory.h"

#include <float.h>

vtkCxxRevisionMacro(vtkcvisHedgeHog, "$LastChangedRevision: 816 $");
vtkStandardNewMacro(vtkcvisHedgeHog);

vtkcvisHedgeHog::vtkcvisHedgeHog()
{
  this->VectorScaleFactor[0] = 0.0;
  this->VectorScaleFactor[1] = 1.0;
  this->ScalingRange[0] = 0.0;
  this->ScalingRange[1] = 1.0;

  this->Crop = 0;
  this->CroppingRange[0] = 0.0;
  this->CroppingRange[1] = 1.0;
}

void vtkcvisHedgeHog::Execute()
{
  vtkDataSet *input= this->GetInput();
  vtkIdType numPts;
  vtkPoints *newPts;
  vtkPointData *pd;
  vtkDataArray *inVectors;
  vtkDataArray *inNormals;
  vtkIdType ptId;
  int i;
  vtkIdType pts[2];
  vtkCellArray *newLines;
  float *x, *v;
  float newX[3];
  vtkPolyData *output = this->GetOutput();
  vtkPointData *outputPD = output->GetPointData();
  vtkDataArray *ptScalars;

  // Initialize
  //
  numPts = input->GetNumberOfPoints();
  pd = input->GetPointData();
  inVectors = pd->GetVectors();
  if ( numPts < 1 )
  {
     vtkErrorMacro(<<"No input data");
     return;
  }
  if ( !inVectors)
  {
     vtkErrorMacro(<<"No vectors in input data");
     return;
  }

  if ( !(ptScalars = pd -> GetScalars()) ) 
  {
     vtkErrorMacro(<<"No scalars in input data");
     return;
  }

  outputPD->CopyAllocate(pd, 2*numPts);
  
  newPts = vtkPoints::New(); newPts->SetNumberOfPoints(2*numPts);
  newLines = vtkCellArray::New();
  newLines->Allocate(newLines->EstimateSize(numPts,2));

  float norm, scale_factor;

  // Loop over all points, creating oriented line
  //
  for (ptId=0; ptId < numPts; ptId++)
  {
     if ( ! (ptId % 10000) ) //abort/progress
     {
	this->UpdateProgress ((float)ptId/numPts);
	if (this->GetAbortExecute())
        {
	   break;
        }
     }
    
    x = input->GetPoint(ptId);
    
    v = inVectors->GetTuple(ptId);

    norm = ptScalars -> GetComponent(ptId, 0);

    if (!Crop || ( (norm > CroppingRange[0]) && norm < CroppingRange[1])) 
    {

       if ( (norm <= FLT_MIN)  || 
	   ((ScalingRange[1] - ScalingRange[0]) <= FLT_MIN) )
       {
	  scale_factor = 0.0;
       }
       else
       {
	  scale_factor = VectorScaleFactor[0] / norm +
	     (VectorScaleFactor[1] - VectorScaleFactor[0]) / 
	     (ScalingRange[1] - ScalingRange[0]);
       }
       
       for (i=0; i<3; i++)
       {
	  newX[i] = x[i] + scale_factor * v[i];
       }
       
       pts[0] = ptId;
       pts[1] = ptId + numPts;;
       
       newPts->SetPoint(pts[0], x);
       newPts->SetPoint(pts[1], newX);
       
       newLines->InsertNextCell(2,pts);
       
       outputPD->CopyData(pd,ptId,pts[0]);
       outputPD->CopyData(pd,ptId,pts[1]);
    }
  } 

  // Update ourselves and release memory
  //
  output->SetPoints(newPts);
  newPts->Delete();

  output->SetLines(newLines);
  newLines->Delete();
}

void vtkcvisHedgeHog::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "VectorScaleFactor: " << this->VectorScaleFactor[0] << " " <<
     this -> VectorScaleFactor[1] << "\n";
  os << indent << "Crop: " << this->Crop << "\n";
  os << indent << "CroppingRange: " << this->CroppingRange[0] << " " <<
     this -> CroppingRange[1] << "\n";
}
