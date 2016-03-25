// Test1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vtkMath.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkOutlineFilter.h>

#include <samraiStructuredPointsReader.h>

int main(int argc, char* argv[])
{
	vtkMath *math;
	math = vtkMath::New();

	vtkStructuredPoints *volume;

#if 0

	volume = vtkStructuredPoints::New();

	int dim[] = {5, 10, 15};
	volume -> SetDimensions(dim);

	int numScalars;

	numScalars = 4*9*14;
    
	vtkScalars *cellScalars;

	cellScalars = vtkScalars::New();

	cellScalars -> SetNumberOfScalars(numScalars);

	for(int i = 0; i < numScalars; i++) 
	{
		cellScalars -> SetScalar(i, math -> Random(0,1));
	}

    volume -> GetCellData() -> SetScalars(cellScalars);
#else

	samraiStructuredPointsReader *samraiVolume;
	samraiVolume = samraiStructuredPointsReader::New();

	samraiVolume -> SetFileName("u:\\users\\ssmith\\SAMRAI\\VTK\\tests\\amr3d.0000.vis");

	samraiVolume -> Execute();

   samraiVolume ->ComputeBounds();
   samraiVolume ->ComputeRange();

   for(int i = 0; i < samraiVolume -> GetNumberOfVariables(); i++)
   {
      cout << ":" << samraiVolume -> GetVariableName(i) << ":\n";
   }

   // samraiVolume -> SetVariableName("Pressure");

   volume = samraiVolume -> GetPatch(0);
#endif

	vtkPlane *plane;
	plane = vtkPlane::New();

	plane -> SetOrigin(volume -> GetCenter());
	
	float normal[] = {0, 1, 1};
	plane -> SetNormal(normal);


	vtkCutter * planeCut;
	planeCut = vtkCutter::New();

	planeCut -> SetInput(volume);
	planeCut -> SetCutFunction(plane);

	vtkPolyDataMapper *mapper;
	mapper = vtkPolyDataMapper::New();

	mapper -> SetInput(planeCut -> GetOutput());
	mapper -> GlobalImmediateModeRenderingOn();

   mapper -> SetScalarRange(samraiVolume -> GetScalarRange());

	vtkActor *cutActor;
	cutActor = vtkActor::New();
	cutActor -> SetMapper(mapper);

	vtkOutlineFilter *outline;
	outline = vtkOutlineFilter::New();
	outline -> SetInput(volume);

	vtkPolyDataMapper *outlineMapper;
	outlineMapper = vtkPolyDataMapper::New();
	outlineMapper -> SetInput(outline -> GetOutput());

	vtkActor *outlineActor;
	outlineActor = vtkActor::New();
	outlineActor -> SetMapper(outlineMapper);

	vtkRenderer *renderer = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(renderer);
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

	renderer -> AddActor(cutActor);
	renderer -> AddActor(outlineActor);
	renderer -> SetBackground(0.1, 0.2, 0.4);
	renWin -> SetSize(300,300);


    // interact with data
	renWin->Render();
	iren->Start();

	return 0;
}
