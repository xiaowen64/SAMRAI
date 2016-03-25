//
// File:        SamraiReader.cpp
// Package:     Vizamrai
// Copyright:   (c) 1997-2000 The Regents of the University of California
// Revision:    $Revision: 45 $
// Modified:    $Date: 2004-12-09 08:25:19 -0800 (Thu, 09 Dec 2004) $
// Description: A utility to read a Vizamrai output file and dump values
//              to the screen
//

#include <stdio.h>
#include <iostream.h>

#include <stdlib.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#define MAX_STRING 2056
#define SAMRAI_NULL_VALUE 1.0E30

#define VTK_LARGE_FLOAT 1.0e+38F
#define VTK_FLOAT_MIN -VTK_LARGE_FLOAT
#define VTK_FLOAT_MAX VTK_LARGE_FLOAT

int xdr_read_string(XDR *stream, char **data, int max_length)
{
   if (!xdr_string(stream, data, max_length))
   {
      return 1;
   }
   return 0;
}

// xdr read routines from Scott 
int xdr_read_int(XDR *stream, int *data, const int n)
{
   if (!xdr_vector(stream, (char *) data, n,
                   sizeof(int), (xdrproc_t) xdr_int))
   {
      return 1;
   }
   return 0;
}

int xdr_read_float(XDR *stream, float *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(float), (xdrproc_t) xdr_float))
   {
      return 1;
   }
   return 0;
}

int xdr_read_double(XDR *stream, double *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(double), (xdrproc_t) xdr_double))
   {
      return 1;
   }
   return 0;
}

int xdr_write_string(XDR *stream, char **data, int max_length)
{
   if (!xdr_string(stream, data, max_length))
   {
      return 1;
   }
   return 0;
}

// xdr read routines from Scott 
int xdr_write_int(XDR *stream, int *data, const int n)
{
   if (!xdr_vector(stream, (char *) data, n,
                   sizeof(int), (xdrproc_t) xdr_int))
   {
      return 1;
   }
   return 0;
}

int xdr_write_float(XDR *stream, float *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(float), (xdrproc_t) xdr_float))
   {
      return 1;
   }
   return 0;
}

int xdr_write_double(XDR *stream, double *data, const int n)
{
   if (!xdr_vector(stream, (char*) data, n,
                   sizeof(double), (xdrproc_t) xdr_double))
   {
      return 1;
   }
   return 0;
}



#if 0
void ComputeRange()
{
   int variable;
   int patch;
   
   for(variable = 0; variable < NumVariables; variable ++)
   {
      ScalarRange[variable][0] = VTK_FLOAT_MAX;
      ScalarRange[variable][1] = VTK_FLOAT_MIN;
   }

   for(patch = 0; patch < NumPatches; patch++)
   {
      for(variable = 0; variable < NumVariables; variable ++)
      {
	 float *patchrange;
	 Scalars[patch][variable] -> ComputeRange();
	 
	 patchrange = Scalars[patch][variable] -> GetRange();

	 if ( ScalarRange[variable][0] > patchrange[0]) 
	    ScalarRange[variable][0] = patchrange[0];

	 if ( ScalarRange[variable][1] < patchrange[1]) 
	    ScalarRange[variable][1] = patchrange[1];
      }
   }

   // If no range then create a small artificial one to prevent errors
   // Could be fixed elsewhere but patching this here seemed to be easier
   for(variable = 0; variable < NumVariables; variable ++)
   {
      if ( ScalarRange[variable][0] == ScalarRange[variable][1] )
	    ScalarRange[variable][1] += 0.000001;
   }
}
#endif

int main(int argc, char **argv)
{
   char * FileName;
   int FileFormat;
   
   int NumPatches;
   
   char *VariableName;
   int Variable;
   
   int NumVariables;
   char **VariableNames;
   
   int *PatchLevel;

   int *Size;
   
   int NumLevels;
   
   int Scales[3];
   
   float **ScalarRange;
   
   float Bounds[6];

   float Spacing[3];
   
   float Scaling;

   unsigned int **StartingPosition;

   int NumDims = 3;

   int DataFormat;

   int DataElementSize;
   
   FILE *file;

   FILE *out_file;

   int dim[3];

   int max_level = 0;

   // Should these be of type float?  or double?  both?

   double delta[3];
   double corner[3];

   XDR stream;
   XDR out_stream;

   int i;
   int ii;

   int patch;

   // Sample Reading 
   if ( (file = fopen(argv[1], "r")) == NULL )
   {
      cerr << "Can't open file" << argv[1] << endl;
      perror("Error reported from fopen");
      return 1;
   }


   if ( (out_file = fopen("test.vis", "wb")) == NULL )
   {
      cerr << "Can't open file" << argv[1] << endl;
      perror("Error reported from fopen");
      return 1;
   }

   xdrstdio_create(&out_stream, out_file, XDR_ENCODE);


   int FileFormat = 5;
   xdr_write_int(&out_stream, &FileFormat, 1);

   cout << "FileFormat: " << FileFormat << endl;

   if (FileFormat > 6) 
   {
      cout << "File format is not supported" << endl;
      exit(1);
   }

   double time;
   time = 0.0;
   xdr_write_double(&out_stream, &time, 1);
   cout << "Plot Time: " << time << endl;

   // Read in the number of boxes
   NumPatches = 1;
   xdr_write_int(&out_stream, &NumPatches, 1);
   cout << "NumPatches: " << NumPatches << endl;

   // Default to Float data
   DataFormat=1;
   xdr_write_int(&out_stream, &DataFormat, 1);

   switch(DataFormat)
   {
      case 0:
      {
	 DataElementSize = 8;
	 cout << "DataFormat: Doubles " << endl;
	 break;
      }
      case 1:
      {
	 DataElementSize = 4;
	 cout << "DataFormat: floats  " << endl;
	 break;
      }
      
   }

   NumVariables = 1;
   xdr_write_int(&out_stream, &NumVariables, 1);
   cout << "NumVariables: " << NumVariables << endl;

   VariableNames = (char **)calloc(NumVariables, sizeof(char *));
   for ( int i = 0; i < NumVariables; i++)
   {
      VariableNames[i] = "Population";
      xdr_write_string(&out_stream, &VariableNames[i], MAX_STRING);
      cout << "VariableNames[" << i << "]: " << VariableNames[i] << endl;
   }

   int *scale;
   NumLevels = 1;
   xdr_write_int(&out_stream, &NumLevels, 1);

   cout << "NumLevels: " << NumLevels << endl;
	 
   scale = (int *)malloc(sizeof(int)*(NumLevels-1)*3);
	 
   xdr_write_int(&out_stream, scale, (NumLevels-1)*3);

   for(ii = 0; ii < NumLevels-1; ii++)
   {
      cout << "scale[" << ii << "]: <"  
	   << scale[ii*3] << ","
	   << scale[ii*3+1] << ","
	   << scale[ii*3+2] << ">" << endl;
   }
	 
   free(scale);

   PatchLevel = (int *)calloc(NumPatches, sizeof(int));

   for( i = 0; i < NumDims; i++) 
   {
      Spacing[i] = VTK_FLOAT_MAX;
   }

   int NumberOfPatchBoundaries;
   int BoundaryLevel;

   NumberOfPatchBoundaries = 1;
   xdr_write_int(&out_stream, &NumberOfPatchBoundaries, 1);
   for(patch = 0; patch < NumberOfPatchBoundaries; patch ++)
   {
      double d_bounds[6];


      fscanf(file, "%d %d %d", &dim[0], &dim[1], &dim[2]);
      fscanf(file, "%f %f %f", &Spacing[0], &Spacing[1], &Spacing[2]);
      
      d_bounds[0] = 0.0;
      d_bounds[1] = 0.0;
      d_bounds[2] = 0.0;

      d_bounds[3] = dim[0] * Spacing[0];
      d_bounds[4] = dim[1] * Spacing[1];
      d_bounds[5] = dim[2] * Spacing[2];

      xdr_write_double(&out_stream, d_bounds, 6);
      cout << "Patch Boundary: "  << patch << endl;
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "    Lower[" << ii << "]: " << d_bounds[ii] << endl;
      }

      for(ii=3; ii < 6; ii++) 
      {
	 cout << "    Upper[" << ii << "]: " << d_bounds[ii] << endl;
      }

      xdr_read_int(&stream, &BoundaryLevel, 1);
      xdr_write_int(&out_stream, &BoundaryLevel, 1);
   }

   for(patch = 0; patch < NumPatches; patch++)
   {
      int lower[3];
      int upper[3];
      int size;

      PatchLevel[patch] = 0;
      xdr_write_int(&out_stream, &PatchLevel[patch], 1);
      cout << "PatchLevel[ " << patch << "]: " << PatchLevel[patch] 
	   << endl;
      max_level = PatchLevel[patch] > max_level ? 
	 PatchLevel[patch] : max_level;
   
      lower[0] = 0;
      lower[1] = 0;
      lower[2] = 0;

      upper[0] = dim[0];
      upper[1] = dim[1];
      upper[2] = dim[2];

      xdr_write_int(&out_stream, lower, 3);
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "Lower[" << ii << "]: " << lower[ii] << endl;
      }

      xdr_write_int(&out_stream, upper, 3);
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "Upper[" << ii << "]: " << upper[ii] << endl;
      }

      delta[0] = Spacing[0];
      delta[1] = Spacing[1];
      delta[2] = Spacing[2];

      xdr_write_double(&out_stream, delta, 3);
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "Delta[" << ii << "]: " << delta[ii] << endl;
      }

      corner[0] = delta[0] / 2.0;
      corner[1] = delta[1] / 2.0;
      corner[2] = delta[2] / 2.0;

      xdr_write_double(&out_stream, corner, 3);
      for(ii=0; ii < 3; ii++) 
      {
	 cout << "Corner[" << ii << "]:" << corner[ii] << endl;
      }

      dim[0] = upper[0]-lower[0]+1;
      dim[1] = upper[1]-lower[1]+1;
      dim[2] = upper[2]-lower[2]+1;


      // Number of Cells
      size = dim[0]*dim[1]*dim[2];

      float ar[3];

      switch (DataFormat) 
      {
	 // Read doubles
	 case 0:
	 {
	    // Don't use doubles
	    double *temp = new double[size];
	    delete[] temp;
	    break;
	 }

	 // Read floats
	 case 1:
	 {
	    float *temp = new float[size*NumVariables];

	    // SGS read in the floats
	    xdr_write_float(&out_stream, temp, size*NumVariables);

	    delete[] temp;
	    break;
	 }
      }

      // End of patch loop
   }

   fflush(out_file);
   xdr_destroy(&out_stream);
   fclose(file);
   fclose(out_file);

   return 0;
}

