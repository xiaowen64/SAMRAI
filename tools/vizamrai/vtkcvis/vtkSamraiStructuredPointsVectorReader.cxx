//
// File:        vtkSamraiStructuredPointsVectorReader.cxx
// Package:     Vizamrai
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision$
// Modified:    $Date$
// Description: File reader for Vizamrai files
//

#include "vtkSamraiStructuredPointsVectorReader.h"

#include <vtkObjectFactory.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#define SCALE_SIZE 1e6


// This is an arbitrary length for the maximum length of an
// XDR string.  Need to have some max length to provide to
// XDR and this seems sufficient.  Famous last words...
#define VARIABLE_NAME_MAX_LENGTH 5000

// This is a value which should not be visualized.  This
// is currently not used.
#define SAMRAI_NULL_VALUE 1.0E30

vtkCxxRevisionMacro(vtkSamraiStructuredPointsVectorReader, "$Revision$");
vtkStandardNewMacro(vtkSamraiStructuredPointsVectorReader);


vtkSamraiStructuredPointsVectorReader::vtkSamraiStructuredPointsVectorReader()
{
   FileName = NULL;

   Size = NULL;
   StartingPosition = NULL;
   
   DataNumPatches = 0;
   DisplayNumPatches = 0;

   Variable = 0; 
   VariableName = NULL;
   VariableNames = NULL;

   Index = 0;

   DataPatchLevel = NULL;
   DisplayPatchLevel = NULL;

   NumLevels = 0;

   Time = 0.0;

   DataStructuredPoints = NULL;
   DisplayStructuredPoints = NULL;

   DisplayScalars = NULL;

   DataDim = NULL;
   DisplayDim = NULL;
   DataSpacing = NULL;

   Refinement = NULL;
   Multiple = NULL;

   File = NULL;

   NumberOfPatchBoundaries = 0;
   PatchBoundaries = NULL;
   PatchBoundaryLevel = NULL;

   Stride[0] = Stride[1] = Stride[2] = 1;

   Collection = vtkStructuredPointsCollection::New();

   Scaling = 1.0;

   Cropping = 0;
   CroppingRegionPlanes[0] = 0;
   CroppingRegionPlanes[1] = 1;
   CroppingRegionPlanes[2] = 0;
   CroppingRegionPlanes[3] = 1;
   CroppingRegionPlanes[4] = 0;
   CroppingRegionPlanes[5] = 1;
}

vtkSamraiStructuredPointsVectorReader::~vtkSamraiStructuredPointsVectorReader()
{
   FreeAllocated();
   
   Collection -> Delete();
}

const char *vtkSamraiStructuredPointsVectorReader::GetClassName() 
{
      return "vtkSamraiStructuredPointsVectorReader";
}

// Get the type of file 
int vtkSamraiStructuredPointsVectorReader::GetFileFormat() 
{
   return FileFormat;
}

char *vtkSamraiStructuredPointsVectorReader::GetCurrentVariableName() 
{
   return VariableName;
}

void vtkSamraiStructuredPointsVectorReader::SetIndex(int index) 
{
   Index = index;
   // Fetch new data
   SetCurrentVariableName(GetCurrentVariableName());
}

int vtkSamraiStructuredPointsVectorReader::GetDepth(int i)
{
   return Depths[i];
}

// Set the name of the scalar data to extract. If not specified, first 
// scalar data encountered is extracted.
void vtkSamraiStructuredPointsVectorReader::SetCurrentVariableName(char *name) 
{
   int i,j,k;

   int patch;

   for(i = 0; i < NumVariables; i++)
   {
      if(strcmp(name, VariableNames[i]) == 0)
         break;
   }

   if(i >= NumVariables)
   {
      // Default to first variable name
      i = 0;
      return;
   }

   if(DisplayScalars) 
   {
      // Remove the old variable data
      for(int patch = 0; patch < DisplayNumPatches; patch++)
      {
	 DisplayScalars[patch] -> Delete();
	 DisplayScalars[patch] = NULL;
      }
   }

   Variable = i;
   VariableName = VariableNames[i];

   DisplayNumPatches = 0;

   for(patch = 0; patch < DataNumPatches; patch++)
   {
      bool inside[3];
      int size = Size[patch];
      float *bounds;

      bounds = DataStructuredPoints[patch] -> GetBounds();

      if( ((bounds[0] >= CroppingRegionPlanes[0]) &&
	   (bounds[0] <= CroppingRegionPlanes[1])) ||
	  ((bounds[1] >= CroppingRegionPlanes[0]) &&
	   (bounds[1] <= CroppingRegionPlanes[1])) )
	 inside[0] = true;
      else
	 inside[0] = false;

      if( ((bounds[2] >= CroppingRegionPlanes[2]) &&
	   (bounds[2] <= CroppingRegionPlanes[3])) ||
	  ((bounds[3] >= CroppingRegionPlanes[2]) &&
	   (bounds[3] <= CroppingRegionPlanes[3])) )
	 inside[1] = true;
      else
	 inside[1] = false;

      if( ((bounds[4] >= CroppingRegionPlanes[4]) &&
	   (bounds[4] <= CroppingRegionPlanes[5])) ||
	  ((bounds[5] >= CroppingRegionPlanes[4]) &&
	   (bounds[5] <= CroppingRegionPlanes[5])) )
	 inside[2] = true;
      else
	 inside[2] = false;

      if( (CroppingRegionPlanes[0] >= bounds[0]) &&
	  (CroppingRegionPlanes[1] <= bounds[1]) )
	 inside[0] = true;

      if( (CroppingRegionPlanes[2] >= bounds[2]) &&
	  (CroppingRegionPlanes[3] <= bounds[3]) )
	 inside[1] = true;

      if( (CroppingRegionPlanes[4] >= bounds[4]) &&
	  (CroppingRegionPlanes[5] <= bounds[5]) )
	 inside[2] = true;
      
      if( !Cropping || (inside[0] && inside[1] && inside[2]))
      {
	 
	 DisplayScalars[DisplayNumPatches] = vtkFloatArray::New();

	 // VTK appears to only allow 3 dim vectors for hedgehogging
	 DisplayScalars[DisplayNumPatches] -> SetNumberOfComponents(3);
	 
	 DisplayScalars[DisplayNumPatches] -> SetNumberOfTuples(size);

	 DisplayStructuredPoints[DisplayNumPatches] = 
	    DataStructuredPoints[patch];
	 
	 vtkFloatArray *data_array = 
	    (vtkFloatArray *)DisplayScalars[DisplayNumPatches];

	 for(int component = 0; component < Depths[Variable]; component++)
	 {
	    xdr_setpos(&Stream, StartingPosition[patch][Variable][component]);

	    int offset = 0;
	 
	    // above here
	    switch (DataFormat) 
	    {
	       // Read doubles
	       case 0:
	       {
		  double *temp = new double[DataDim[patch][0]];
		  
		  for(k = 0; k < DataDim[patch][2]; k += Stride[2]) 
		  {
		     for(j = 0; j < DataDim[patch][1]; j += Stride[1]) 
		     {
			xdr_read_double(&Stream, temp, DataDim[patch][0]);
			for(i = 0; i < DataDim[patch][0]; i += Stride[0]) 
			{
			   data_array -> SetComponent(offset++, 
						      component, temp[i]);
			}
			xdr_setpos(&Stream, xdr_getpos(&Stream) + 
				   (Stride[1]-1) * DataDim[patch][0] 
				   * DataElementSize);
		     }
		     int backup;
		     if ( backup = DataDim[patch][1] % Stride[1] ) 
		     {
			backup = backup - Stride[1];
		     }
		     
		     xdr_setpos(&Stream, xdr_getpos(&Stream) + 
				(DataDim[patch][1] * DataDim[patch][0] * (Stride[2]-1) -
				 backup*DataDim[patch][0])
				* DataElementSize);
		  }
	       
		  delete[] temp;
		  break;
	       }
	    
	       // Read floats
	       case 1:
	       {
		  float *temp = new float[DataDim[patch][0]];
		  
		  for(k = 0; k < DataDim[patch][2]; k += Stride[2]) 
		  {
		     for(j = 0; j < DataDim[patch][1]; j += Stride[1]) 
		     {
			xdr_read_float(&Stream, temp, DataDim[patch][0]);
			for(i = 0; i < DataDim[patch][0]; i += Stride[0]) 
			{
			   
			   data_array -> SetComponent(offset++, 
						      component, temp[i]);
			}
			
			// Check if we are at the end of the data to see how far
			// we should skip.  If at the end then we may not want to 
			// go a full stride, just need to skip over the remainder
			// on this 2D plane.
			if ( (j+Stride[1]) < DataDim[patch][1])
			{
			   xdr_setpos(&Stream, xdr_getpos(&Stream) + 
				      (Stride[1]-1) * DataDim[patch][0] 
				      * DataElementSize);
			}
			else
			{
			   xdr_setpos(&Stream, xdr_getpos(&Stream) + 
				      ((DataDim[patch][1]  - 1) % Stride[1]) * DataDim[patch][0] 
				      * DataElementSize);
			}
		     }
		     xdr_setpos(&Stream, xdr_getpos(&Stream) + 
				(DataDim[patch][1] * DataDim[patch][0] * (Stride[2]-1))
				* DataElementSize);
		  }

		  delete[] temp;
		  break;
	       }
	    }
	 }

	 if (Depths[Variable] == 2) {
	    for(int i = 0; i < size; i++) 
	    {
	       data_array -> SetComponent(i, 2, 0.0);

	    }
	 }
	 
	 DisplayScalars[DisplayNumPatches] -> Modified();

	 // ************************************************************************
	 // FIXME for VTK 4.0!!!!!!!!!!!!!!!!!!!!!!!!!!

	 DisplayStructuredPoints[DisplayNumPatches] -> GetCellData() -> 
	    SetScalars(DisplayScalars[DisplayNumPatches]);
	 
	 DisplayStructuredPoints[DisplayNumPatches] -> GetPointData() -> Update();
	 
	 // Increment number of visible patches
	 DisplayNumPatches++;
      }
   }

   ComputeRange();
}

char **vtkSamraiStructuredPointsVectorReader::GetVariableNames()
{
   return VariableNames;
}

char *vtkSamraiStructuredPointsVectorReader::GetVariableName(int i)
{
   return VariableNames[i];
}

int vtkSamraiStructuredPointsVectorReader::GetNumberOfVariables() 
{
   return NumVariables;
}


int vtkSamraiStructuredPointsVectorReader::GetNumberOfPatches()
{
   return DisplayNumPatches;
}

vtkStructuredPoints *vtkSamraiStructuredPointsVectorReader::GetPatch(int patch)

{
   return DisplayStructuredPoints[patch];
}

int vtkSamraiStructuredPointsVectorReader::GetNumberOfRefinementLevels()
{
   return RefinementLevels;
}

int *vtkSamraiStructuredPointsVectorReader::GetRefinement(int level)
{
   return &Refinement[level*3];
}

int vtkSamraiStructuredPointsVectorReader::GetPatchLevel(int patch)
{
   return DataPatchLevel[patch];
}

int vtkSamraiStructuredPointsVectorReader::GetNumberOfLevels()
{
   return NumLevels;
}

float *vtkSamraiStructuredPointsVectorReader::GetBounds() 
{
   return Bounds;
}

int *vtkSamraiStructuredPointsVectorReader::GetIndexBounds() 
{
   return IndexBounds;
}

float *vtkSamraiStructuredPointsVectorReader::GetIndexScaling() 
{
   return IndexScale;
}


void   vtkSamraiStructuredPointsVectorReader::ComputeBounds()
{
   Bounds[0] = Bounds[2] = Bounds[4] = VTK_FLOAT_MAX;
   Bounds[1] = Bounds[3] = Bounds[5] = VTK_FLOAT_MIN;

   for(int patch = 0; patch < DataNumPatches; patch++)
   {
      float *patchbounds;

      DataStructuredPoints[patch] -> ComputeBounds();
      patchbounds = DataStructuredPoints[patch] -> GetBounds();
      
      switch (FileFormat)
      {
	 case 1:
	 case 2:
	 case 3:
         {
	    PatchBoundaries[patch] -> SetBounds(patchbounds);
	    PatchBoundaryLevel[patch] = DataPatchLevel[patch];
	    break;
	 }
      }

      // Set the mins
      if(Bounds[0] > patchbounds[0]) 
	 Bounds[0] = patchbounds[0];

      if(Bounds[2] > patchbounds[2]) 
	 Bounds[2] = patchbounds[2];

      if(Bounds[4] > patchbounds[4]) 
	 Bounds[4] = patchbounds[4];

      // Set the maxs
      if(Bounds[1] < patchbounds[1]) 
	 Bounds[1] = patchbounds[1];

      if(Bounds[3] < patchbounds[3]) 
	 Bounds[3] = patchbounds[3];

      if(Bounds[5] < patchbounds[5]) 
	 Bounds[5] = patchbounds[5];
   }
}

float *vtkSamraiStructuredPointsVectorReader::GetScalarRange()
{
   return ScalarRange;
}

void vtkSamraiStructuredPointsVectorReader::ComputeRange()
{
   int patch;

   ScalarRange[0] = VTK_FLOAT_MAX;
   ScalarRange[1] = VTK_FLOAT_MIN;

   for(patch = 0; patch < DisplayNumPatches; patch++)
   {
      float *patchrange;
      DisplayScalars[patch] -> ComputeRange(0);
      
      patchrange = DisplayScalars[patch] -> GetRange();
      
      if ( ScalarRange[0] > patchrange[0]) 
	 ScalarRange[0] = patchrange[0];
      
      if ( ScalarRange[1] < patchrange[1]) 
	 ScalarRange[1] = patchrange[1];
   }

   // If no range then create a small artificial one to prevent errors
   // Could be fixed elsewhere but patching this here seemed to be easier
   if ( ScalarRange[0] == ScalarRange[1]) 
   {
      if ( ScalarRange[0] == 0.0 ) 
	 ScalarRange[0] = -1;
      else
	 ScalarRange[0] -= ScalarRange[0]*FLT_EPSILON;
   }
}

float *vtkSamraiStructuredPointsVectorReader::GetSpacing()
{
   return Spacing;
}

float *vtkSamraiStructuredPointsVectorReader::GetCenter()
{
   float *bounds;

   bounds = GetBounds();

   Center[0] = (bounds[1] - bounds[0])/2 + bounds[0];
   Center[1] = (bounds[3] - bounds[2])/2 + bounds[2];
   Center[2] = (bounds[5] - bounds[4])/2 + bounds[4];

   return Center;
}

void vtkSamraiStructuredPointsVectorReader::Scale()
{
   
   ComputeRange();

   for(int patch = 0; patch < DisplayNumPatches; patch++)
   {
      int size = DisplayScalars[patch] -> GetNumberOfTuples();
      
      if( ScalarRange[1] - ScalarRange[0] > 0.0)
      {
	 for(int i = 0; i < size; i++)
	 {
	    float value;
	    // Slow way to do this; should get the array out so don't need
	    // Overhead of function calls
	    value = DisplayScalars[patch] -> GetTuple1(i);
	    
	    value = (value - ScalarRange[0]) / 
	       (ScalarRange[1] - ScalarRange[0]);
	    
	    DisplayScalars[patch] -> SetTuple1(i, value);
	 }
      }
      else
      {
	 // If no range then just set to zero
	 for(int i = 0; i < size; i++)
	 {
	    DisplayScalars[patch] -> SetTuple1(i, 0.0);
	 }
      }
      
      DisplayScalars[patch] -> Modified();
   }
      
   // This is a slow way to do this
   ComputeRange();
}

void vtkSamraiStructuredPointsVectorReader::Execute()
{

   int NumDims = 3;
   
   int max_level = 0;

   // Should these be of type float?  or double?  both?

   double delta[3];
   double corner[3];

   int i;

   int patch;
   int variable;

   // First delete existing data if any
   FreeAllocated();

   for(i = 0; i < NumDims; i++)
   {
      IndexBounds[i*2+0] = VTK_INT_MAX;
      IndexBounds[i*2+1] = VTK_INT_MIN;
   }
   
   // Sample Reading 
   if ( (File = fopen(FileName, "rb")) == NULL )
   {
      vtkErrorMacro(<< "Can't open file" << FileName << "\n");
      return;
   }

   xdrstdio_create(&Stream, File, XDR_DECODE);

   // Read in the first int in file; this can be the number of boxes
   // or the file type (if < 0)
   xdr_read_int(&Stream, &DataNumPatches, 1);

   

   if ( DataNumPatches < 0 ) 
   {
      // Dealing with new format 
      FileFormat = abs(DataNumPatches);
   }
   else
      FileFormat = 0;

   if (FileFormat > 6) 
   {
      vtkErrorMacro(<< "File format for " << FileName << " is unsupported\n");
      return;
   }

   switch (FileFormat) 
   {
      case 5:
      case 6:
      {
	 // Get the time for this data 
	 xdr_read_double(&Stream, &Time, 1);
	 break;
      }
   }

   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 // Read in the number of boxes
	 xdr_read_int(&Stream, &DataNumPatches, 1);
	 break;
      }
   }


   // Default to double for older file formats
   DataFormat=0;
   switch (FileFormat) 
   {
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 // Read storage format for data values
	 xdr_read_int(&Stream, &DataFormat, 1);
	 break;
      }
   }

   
   switch(DataFormat)
   {
      case 0:
      {
	 DataElementSize = 8;
	 break;
      }
      case 1:
      {
	 DataElementSize = 4;
	 break;
      }
      
   }

   xdr_read_int(&Stream, &NumVariables, 1);

   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
         VariableNames = (char **)calloc(NumVariables, sizeof(char *));
         Depths = (int *)calloc(NumVariables, sizeof(int));
         for ( int i = 0; i < NumVariables; i++)
         {
            // Read in the name for the component
            xdr_read_string(&Stream, &VariableNames[i], VARIABLE_NAME_MAX_LENGTH);

	    switch (FileFormat) 
	    {
	       case 1:
	       case 2:
	       case 3:
	       case 4:
	       case 5:
	       {
		  Depths[i] = 1;
		  break;
	       }
	       case 6:
	       {
		  xdr_read_int(&Stream, &Depths[i], 1);
		  break;
	       }
	       break;
	    }
	 }
      }
   }


   switch (FileFormat) {
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 xdr_read_int(&Stream, &RefinementLevels, 1);
	 
	 Refinement = (int *)malloc(sizeof(int)*RefinementLevels*3);
	 Multiple   = (int *)malloc(sizeof(int)*RefinementLevels*3);

	 // Level 0 is always 1
	 Refinement[0] = Refinement[1] = Refinement[2] = 1;
	 
	 xdr_read_int(&Stream, &Refinement[3], (RefinementLevels-1)*3);

	 for(i=0; i < RefinementLevels; i++ )
	 {
	    Multiple[i*3] = Multiple[i*3+1] = Multiple[i*3+2] = 1;
	 }

	 i = RefinementLevels-1;
	 while(i > -1)
	 {
	    int j;
	    j = i - 1;
	    while(j > -1) 
	    {
	       Multiple[j*3] = Multiple[j*3] * Refinement[i*3];
	       Multiple[j*3+1] = Multiple[j*3+1] * Refinement[i*3+1];
	       Multiple[j*3+2] = Multiple[j*3+2] * Refinement[i*3+2];
	       j--;
	    }
	    i--;
	 }
	 break;
      }
   }

   DataPatchLevel = (int *)calloc(DataNumPatches, sizeof(int));
   DisplayPatchLevel = (int *)calloc(DataNumPatches, sizeof(int));

   Size       = (int *)calloc(DataNumPatches, sizeof(int));

   StartingPosition = (unsigned int ***)calloc(DataNumPatches,
					      sizeof(unsigned int **));

   DataStructuredPoints = (vtkStructuredPoints **)
      calloc(DataNumPatches, sizeof(vtkStructuredPoints *));

   DisplayStructuredPoints = (vtkStructuredPoints **)
      calloc(DataNumPatches, sizeof(vtkStructuredPoints *));

   DisplayScalars = (vtkDataArray **)
      calloc(DataNumPatches, sizeof(vtkDataArray *));

   for( i = 0; i < NumDims; i++) 
   {
      Spacing[i] = VTK_FLOAT_MAX;
   }


   switch (FileFormat)
   {
      case 1:
      case 2:
      case 3:
      {
	 NumberOfPatchBoundaries = DataNumPatches;
	 PatchBoundaries = (vtkOutlineSource **)calloc(NumberOfPatchBoundaries,
					       sizeof(vtkOutlineSource *));
	 PatchBoundaryLevel = (int *)calloc(NumberOfPatchBoundaries,
					    sizeof(int));
	 break;
      }
      case 4:
      case 5:
      case 6:
      {
	 xdr_read_int(&Stream, &NumberOfPatchBoundaries, 1);
	 PatchBoundaries = (vtkOutlineSource **)calloc(NumberOfPatchBoundaries,
					       sizeof(vtkOutlineSource *));
	 PatchBoundaryLevel = (int *)calloc(NumberOfPatchBoundaries,
					    sizeof(int));
	 for(patch = 0; patch < NumberOfPatchBoundaries; patch ++)
	 {
	    double d_bounds[6];
	    float  f_bounds[6];
	    
	    PatchBoundaries[patch] = vtkOutlineSource::New();
	    xdr_read_double(&Stream, d_bounds, 6);
	    f_bounds[0] = d_bounds[0];
	    f_bounds[1] = d_bounds[3];

	    f_bounds[2] = d_bounds[1];
	    f_bounds[3] = d_bounds[4];

	    f_bounds[4] = d_bounds[2];
	    f_bounds[5] = d_bounds[5];
	    
	    PatchBoundaries[patch] -> SetBounds(f_bounds);
	    xdr_read_int(&Stream, &PatchBoundaryLevel[patch], 1);
	 }
	 break;
      }
   }

   DataDim = new int *[DataNumPatches];
   DisplayDim = new int *[DataNumPatches];
   DataSpacing = new float *[DataNumPatches];
   
   for(patch = 0; patch < DataNumPatches; patch++)
   {
      int lower[3];
      int upper[3];
      int datasize;
      
      DataStructuredPoints[patch] = vtkStructuredPoints::New();

      StartingPosition[patch] = (unsigned int **)calloc(NumVariables, 
						  sizeof(unsigned int*));
      switch (FileFormat) 
      {
         case 1:
	 case 2:
	 case 3:
         case 4:
         case 5:
         case 6:
         {
            xdr_read_int(&Stream, &DataPatchLevel[patch], 1);
            max_level = DataPatchLevel[patch] > max_level ? 
               DataPatchLevel[patch] : max_level;
            break;
         }
      }
      
      
      xdr_read_int(&Stream, lower, 3);
      xdr_read_int(&Stream, upper, 3);


      for(i = 0; i < NumDims; i++)
      {
	 int background;
	 background = lower[i] * Multiple[DataPatchLevel[patch]*3+i];
	 if (background < IndexBounds[i*2+0]) {
	    IndexBounds[i*2+0] = background;
	 }

	 // Need to adjust so we stop at the upper boundary
	 // |++++|
	 //     ^
	 // here
	 // |++++|
	 //  ^
	 // not here
	 background = (upper[i]+1) * Multiple[DataPatchLevel[patch]*3+i] - 1;
	 if (background > IndexBounds[i*2+1]) {
	    IndexBounds[i*2+1] = background;
	 }
      }
      
      switch (FileFormat) 
      {
         case 1:
         case 2:
	 case 3:
         case 4:
         case 5:
         case 6:
         {
            xdr_read_double(&Stream, delta, 3);
            xdr_read_double(&Stream, corner, 3);
            break;
         }
      }
      
      DataDim[patch] = new int[3];
      DataDim[patch][0] = upper[0]-lower[0]+1;
      DataDim[patch][1] = upper[1]-lower[1]+1;
      DataDim[patch][2] = upper[2]-lower[2]+1;

      DisplayDim[patch] = new int[3];
      DisplayDim[patch][0] = DataDim[patch][0] / Stride[0] + 
	 ( (DataDim[patch][0] % Stride[0]) > 0) ;
      DisplayDim[patch][1] = DataDim[patch][1] / Stride[1] +
	 ( (DataDim[patch][1] % Stride[1]) > 0) ;
      DisplayDim[patch][2] = DataDim[patch][2] / Stride[2] +
	 ( (DataDim[patch][2] % Stride[2]) > 0) ;
      
      // Number of Cells
      Size[patch] = DisplayDim[patch][0]*DisplayDim[patch][1]*
	 DisplayDim[patch][2];

      datasize = DataDim[patch][0]*DataDim[patch][1]*
	 DataDim[patch][2];
      
      DataStructuredPoints[patch] -> SetDimensions(DisplayDim[patch][0]+1, 
					       DisplayDim[patch][1]+1, 
					       DisplayDim[patch][2]+1);
      
      DataSpacing[patch] = new float[3];
      float DisplaySpacing[3];
      
      // Find the smallest delta spacing in all patches, use that for
      // spacing for the entire dataset
      for( i = 0; i < NumDims; i++) 
      {
         DataSpacing[patch][i] = delta[i];
	 DisplaySpacing[i] = (delta[i] * DataDim[patch][i]) 
	    / DisplayDim[patch][i];
	 
	 if ( DataSpacing[patch][i] < Spacing[i] ) 
	    Spacing[i] = DataSpacing[patch][i];
      }

      DataStructuredPoints[patch] -> SetSpacing(DisplaySpacing);
      
      float origin[3];
      
      for( i = 0; i < NumDims; i++) 
         origin[i] = lower[i] * delta[i] + corner[i] - delta[i] * 0.5;

      DataStructuredPoints[patch] -> SetOrigin(origin);

      switch (FileFormat)
      {
	 case 1:
	 case 2:
	 case 3:
         {
	    PatchBoundaries[patch] = vtkOutlineSource::New();
	    break;
	 }
      }

      int num_items = 0;
      for(variable = 0; variable < NumVariables; variable++)
      {
	 StartingPosition[patch][variable] = (unsigned int *)calloc(Depths[variable], 
						  sizeof(unsigned int));
	 for(i = 0; i < Depths[variable]; i++) 
	 {
	    StartingPosition[patch][variable][i] = xdr_getpos(&Stream) + num_items*datasize*DataElementSize;
	    num_items++;
	 }
      }

      // Skip over data to next patch
      xdr_setpos(&Stream, xdr_getpos(&Stream) + num_items * datasize * DataElementSize);
   } // for patch

   NumLevels = max_level + 1;

   // Compute the bounding box for the unscaled domain
   // SGS performance can be improved by doing this calc manually
   ComputeBounds();

   // Compute the conversion factor between background grid index
   // space values and real world coordinates.
   {
      float *patch0_spacing;
      patch0_spacing = DataSpacing[0];
      
      int level;
      level = GetPatchLevel(0);

      for (i = 0; i < NumDims; i++)
      {
	 IndexScale[i] = (patch0_spacing[i] * Multiple[(RefinementLevels-1)*3+i]) /
	    Multiple[(level)*3+i];
      }
   }

   // Scale the coordinates so numbers can be handled by VTK.  
   // VTK was having difficulties with really small numbers

   // Determine scale factor based on bounds.  This should
   // be constant across a given problem run so camera 
   // angles and such should work as users scale the 
   // number of cells.
   Scaling = VTK_FLOAT_MIN;
   for(i = 0; i < NumDims; i++)
   {
      if ((Bounds[i*2+1] - Bounds[i*2]) > Scaling)
	 Scaling = (Bounds[i*2+1] - Bounds[i*2]);
   }
   Scaling = SCALE_SIZE /Scaling;


   // Apply Scaling factor to coordinate values
   for(patch = 0; patch < DataNumPatches; patch++)
   {
      float *origin;
      float *spacing;
      
      float new_origin[3];
      float new_spacing[3];
      
      origin = DataStructuredPoints[patch] -> GetOrigin();
      spacing = DataStructuredPoints[patch] -> GetSpacing();

      for(i = 0; i < NumDims; i++)
      {
	 new_origin[i] = origin[i] * Scaling;
	 new_spacing[i] = spacing[i] * Scaling;
      }

      DataStructuredPoints[patch] -> SetOrigin(new_origin);
      DataStructuredPoints[patch] -> SetSpacing(new_spacing);
   }

   for( i = 0; i < NumDims; i++) 
   {
      Spacing[i] = Spacing[i] * Scaling;
   }


   switch (FileFormat)
   {
      case 4:
      case 5:
      case 6:
	 for(patch = 0; patch < NumberOfPatchBoundaries; patch ++)
	 {
	    float  *f_bounds;
	    f_bounds = PatchBoundaries[patch] -> GetBounds();
	    for(i = 0; i < NumDims*2; i++) 
	    {
	       f_bounds[i] *= Scaling;
	    }
	    PatchBoundaries[patch] -> SetBounds(f_bounds);
	 }
	 break;
   }

   // Compute the bounding box for the domain with the new 
   // scaling
   ComputeBounds();

   // SGS Want to add search for oldname here so we don't 
   // Load in something we don't want ?
   // SetCurrentVariableName(VariableNames[0]); 
}

void vtkSamraiStructuredPointsVectorReader::FreeAllocated()
{
   if(Size) 
   {
      free(Size);
      Size = NULL;
   }

   if(StartingPosition) 
   {
      for(int patch = 0; patch < DataNumPatches; patch++)
      {
	 for(int variable = 0; variable < NumVariables; variable++)
	 {
	    free(StartingPosition[patch][variable]);
	 }
	 free(StartingPosition[patch]);
      }
      free(StartingPosition);
      StartingPosition = NULL;
   }

   if(DataPatchLevel) 
   {
      free(DataPatchLevel);
      DataPatchLevel = NULL;
   }

   if(DisplayPatchLevel) 
   {
      free(DisplayPatchLevel);
      DisplayPatchLevel = NULL;
   }

   if(Refinement)
   {
      free(Refinement);
      Refinement = NULL;
   }

   if(Multiple)
   {
      free(Multiple);
      Multiple = NULL;
   }

   if(DisplayScalars) 
   {
      for(int patch = 0; patch < DisplayNumPatches; patch++)
      {
	 DisplayScalars[patch]-> Delete();
	 DisplayScalars[patch] = NULL;
      }
      DisplayNumPatches = 0;

      free(DisplayScalars);
      DisplayScalars = NULL;
   }

   for(int patch = 0; patch < DataNumPatches; patch++)
   {

      if(DataDim[patch])
	 delete[] DataDim[patch];
      if(DisplayDim[patch])
	 delete[] DisplayDim[patch];

      if(DataSpacing[patch])
	 delete[] DataSpacing[patch];

      DataStructuredPoints[patch] -> Delete();
   }

   if(DataDim) 
   {
      delete[] DataDim;
      DataDim = NULL;
   }

   if(DisplayDim) 
   {
      delete[] DisplayDim;
      DisplayDim = NULL;
   }

   if(DataSpacing) 
   {
      delete[] DataSpacing;
      DataSpacing = NULL;
   }
   

   if (DataStructuredPoints)
   {
      free(DataStructuredPoints);
      DataStructuredPoints = NULL;
   }

   if (DisplayStructuredPoints)
   {
      free(DisplayStructuredPoints);
      DisplayStructuredPoints = NULL;
   }
 
   switch (FileFormat) 
   {
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      {
	 if (VariableNames) 
	 {
	    free(VariableNames);
	    VariableNames = NULL;
	 }
	 break;
      }
   }

   if(PatchBoundaries) 
   {
      for(int patchboundary = 0; patchboundary < NumberOfPatchBoundaries; 
	  patchboundary++)
      {
	 if (PatchBoundaries[patchboundary])
	 {
	    PatchBoundaries[patchboundary] -> Delete();
	 }
      }

      free(PatchBoundaries);
      PatchBoundaries = NULL;
      free(PatchBoundaryLevel);
      PatchBoundaryLevel = NULL;
   }

   if(File) 
   {
      xdr_destroy(&Stream);
      fclose(File);
      File = NULL;
   }
}


// Return the current collection.  This collection should probably be 
// constructed somewhere else and just returned here to improve efficiency.
vtkStructuredPointsCollection *vtkSamraiStructuredPointsVectorReader::GetCollection()
{
   // SGS this is kinda stupid why are we rebuilding the list each 
   // time and updating it?
   Collection -> RemoveAllItems();

   for(int patch = 0; patch < DisplayNumPatches; patch++)
   {
      DisplayStructuredPoints[patch] -> Update();
      Collection -> AddItem(DisplayStructuredPoints[patch]);
   }

   return Collection;
}
   
int vtkSamraiStructuredPointsVectorReader::GetNumberOfPatchBoundaries() 
{
   return NumberOfPatchBoundaries;
}

vtkOutlineSource *vtkSamraiStructuredPointsVectorReader::GetPatchBoundary(int patch)
{
   return PatchBoundaries[patch];
}

int vtkSamraiStructuredPointsVectorReader::GetPatchBoundaryLevel(int patch)
{
   return PatchBoundaryLevel[patch];
}
