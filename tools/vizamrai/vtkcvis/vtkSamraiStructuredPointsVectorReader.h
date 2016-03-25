//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/tools/vizamrai/vtkcvis/vtkSamraiStructuredPointsVectorReader.h $
// Package:     Vizamrai
// Copyright:   (c) 1997-2000 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: File reader for Vizamrai files
//


#ifndef __vtkSamraiStructuredPointsVectorReader_h
#define __vtkSamraiStructuredPointsVectorReader_h

#include <vtkSystemIncludes.h>

#ifdef VTK_USE_ANSI_STDLIB
#include <iomanip>
#else
#include <iomanip.h>
#endif

#ifdef _MSC_VER
#include <winsock.h>
// This is needed by the XDR routines
#define DllExport __declspec( dllimport )
#endif

#include <rpc/types.h>
#include <rpc/xdr.h>
#include <float.h>

#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsCollection.h>
#include <vtkFloatArray.h>
#include <vtkOutlineSource.h>

int xdr_read_string(XDR *stream, char **data, int max_length);
int xdr_read_int(XDR *stream, int *data, const int n);
int xdr_read_float(XDR *stream, float *data, const int n);
int xdr_read_double(XDR *stream, double *data, const int n);

class VTK_EXPORT vtkSamraiStructuredPointsVectorReader : public vtkObject
{
public:

  vtkTypeRevisionMacro(vtkSamraiStructuredPointsVectorReader,vtkObject);

  static vtkSamraiStructuredPointsVectorReader *New();
   
  const char *GetClassName();
   
   // Description:
   // Specify file name of samari data file to read.
   vtkGetStringMacro(FileName);
   vtkSetStringMacro(FileName);

   vtkGetMacro(Time, double);
   vtkSetMacro(Time, double);

   vtkGetMacro(Scaling, float);
   vtkSetMacro(Scaling, float);

   vtkGetVector3Macro(Stride, int);
   vtkSetVector3Macro(Stride, int);

   // Description:
   // Turn On/Off orthogonal cropping. (Clipping planes are
   // perpendicular to the coordinate axes.)
   vtkSetMacro(Cropping,int);
   vtkGetMacro(Cropping,int);
   vtkBooleanMacro(Cropping,int);

   // Description:
   // Set/Get the Cropping Region Planes ( xmin, xmax, ymin, ymax, zmin, zmax )
   vtkSetVector6Macro( CroppingRegionPlanes, float );
   vtkGetVectorMacro(  CroppingRegionPlanes, float, 6 );
   
   // Description:
   // Get the type of file (VTK_ASCII or VTK_BINARY).
   int GetFileFormat();
   
   // Description:
   // Set the name of the scalar data to extract. If not specified, first 
   // scalar data encountered is extracted.
   char *GetCurrentVariableName();
   void  SetCurrentVariableName(char *name);
   
   char **GetVariableNames();
   char * GetVariableName(int i);
   int    GetNumberOfVariables();

   vtkGetMacro(Index, int);

   void SetIndex(int index);

   int GetDepth(int i);
   
   int                  GetNumberOfPatches();
   vtkStructuredPoints *GetPatch(int patch);

   int GetPatchLevel(int patch);
   int GetNumberOfLevels();

   int *GetRefinement(int level);
   int  GetNumberOfRefinementLevels();
   
   float *GetBounds();

   void   ComputeBounds();

   int *GetIndexBounds();

   float *GetIndexScaling();
   
   float *GetScalarRange();
   void ComputeRange();
   
   float *GetSpacing();

   float *GetCenter();

   void Execute();

   void Scale();

   vtkStructuredPointsCollection *GetCollection();

   int GetNumberOfPatchBoundaries();

   vtkOutlineSource *GetPatchBoundary(int patch);

   int GetPatchBoundaryLevel(int patch);
   
protected:

   vtkSamraiStructuredPointsVectorReader();
   ~vtkSamraiStructuredPointsVectorReader();

   char * FileName;
   int FileFormat;
   
   int DataNumPatches;
   int DisplayNumPatches;
   
   // The current working variable name and index
   char *VariableName;
   int Variable;
   int Index;
   
   int NumVariables;
   char **VariableNames;

   int *Depths;
   
   int *DataPatchLevel;
   int *DisplayPatchLevel;

   int *Size;

   int   Cropping;
   float CroppingRegionPlanes[6];
   
   int NumLevels;
   
   int Scales[3];

   int DataFormat;
   int DataElementSize;
   
   float ScalarRange[2];
   float Bounds[6];

   float Center[3];

   float Spacing[3];

   float IndexScale[3];
   int IndexBounds[6];
   
   double Time;

   float Scaling;

   int Stride[3];
   int **DataDim;
   int **DisplayDim;
   
   float **DataSpacing;

   int *Refinement;
   int *Multiple;
   int RefinementLevels;

   unsigned int ***StartingPosition;

   XDR Stream;
   
   FILE *File;
   
   vtkStructuredPoints **DataStructuredPoints;
   vtkStructuredPoints **DisplayStructuredPoints;

   vtkDataArray **DisplayScalars;

   vtkStructuredPointsCollection *Collection;

   int NumberOfPatchBoundaries;
   vtkOutlineSource **PatchBoundaries;
   int *PatchBoundaryLevel;

   void FreeAllocated();
};

#endif



