//
// File:        EmbeddedBoundaryShapePolygonX.C
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 609 $
// Modified:    $Date: 2005-09-13 15:15:49 -0700 (Tue, 13 Sep 2005) $
// Description: Polygon embedded boundary shape
//              
// 

#ifndef included_EmbeddedBoundaryShapePolygon_C
#define included_EmbeddedBoundaryShapePolygon_C

#include "EmbeddedBoundaryShapePolygon.h"

#include "tbox/IEEE.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

#ifdef DEBUG_NO_INLINE
#include "EmbeddedBoundaryShapePolygon.I"
#endif

namespace SAMRAI {
   namespace appu {

template<int DIM>
EmbeddedBoundaryShapePolygon<DIM>::EmbeddedBoundaryShapePolygon(
   const string& object_name,
   tbox::Pointer<tbox::Database> input_db)
{
   d_object_name = object_name;

   d_height = tbox::IEEE::getSignalingNaN();
   d_eps    = tbox::IEEE::getDBL_EPSILON();
   
   getFromInput(input_db);
}

template<int DIM>      
EmbeddedBoundaryShapePolygon<DIM>::~EmbeddedBoundaryShapePolygon<DIM>()
{  
}


template<int DIM>
void EmbeddedBoundaryShapePolygon<DIM>::printClassData(ostream& os) const
{

   os << endl;
   os << "d_object_name = " << d_object_name << endl;
   

   os << "d_eps = " << d_eps << endl;
   os << "d_num_vertices = " << d_num_vertices << endl;
   
   for (int j = 0; j < d_num_vertices; j++) {

      os << "d_vx[" << j << "] = " << d_vx[j] << "\t"
         << "d_vy[" << j << "] = " << d_vy[j] << endl;
   }
   
   if (DIM == 3) {
      os << "height (z) = " << d_height << endl;
   }

   os << endl;
}

template<int DIM>
void EmbeddedBoundaryShapePolygon<DIM>::getFromInput(
   tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   /*
    * Read in coordinates of the nodes of the polygon in X,Y space.
    */
   tbox::Pointer<tbox::Database> vertices_db = db->getDatabase("vertices");
   tbox::Array<string> vertices = vertices_db->getAllKeys();

   d_num_vertices = vertices.getSize();

   d_vx.resizeArray(d_num_vertices);
   d_vy.resizeArray(d_num_vertices);

   int i = 0;

   if (d_num_vertices < 3) {
      TBOX_ERROR(d_object_name << ": " 
                 << "\n You must supply at least 3 vertices to define "
                 << "a polygon." << endl);
   }
   
   tbox::Array<double> vertices_temp;
   for (i=0 ; i < d_num_vertices; i++) {
      string name = vertices[i];

      vertices_temp = vertices_db->getDoubleArray(name);

      if (vertices_temp.getSize() < 2) {
         TBOX_ERROR(d_object_name << ": "
                    << "\ninsufficient entries for 'vertices[ "
                    << i << "]'"
                    << "\n   required size = " << DIM
                    << "\n   supplied size = " << vertices_temp.getSize()
                    << endl);
      }

      d_vx[i] = vertices_temp[0];
      d_vy[i] = vertices_temp[1];
   }


   // We check to make sure that the poly is convex. 
   int ip1, ip2;
   bool counter_clockwise = false;
   bool clockwise = false;
   for (i=0; i < d_num_vertices; i++) {
      
      double v1[3], v2[3], cp[3];
      
      ip1 = (i+1) % d_num_vertices;
      ip2 = (i+2) % d_num_vertices;

      v1[0] = d_vx[ip1] - d_vx[i];
      v1[1] = d_vy[ip1] - d_vy[i];
      v1[2] = 0.;
      v2[0] = d_vx[ip2] - d_vx[ip1];
      v2[1] = d_vy[ip2] - d_vy[ip1];
      v2[2] = 0.;

      crossProduct(cp, v1, v2);
      if (cp[2] > 0)
         counter_clockwise = true;
      else
         clockwise = true;

      if (counter_clockwise && clockwise) // Polygon must not be convex
         TBOX_ERROR(d_object_name << ": "
                    << "\nPolygon must be convex"
                    << "correct vertices in input file."
                    << endl);
            
   }
   
      
   if (DIM == 3) {
      /*
       * MUST supply a height in 3D.
       */
      if (db->keyExists("height")) {
         d_height = db->getDouble("height");
      } else {
         TBOX_ERROR(d_object_name << ": " 
                    << "\n'height' entry not supplied - in 3D you must supply "
                    << "a height." << endl);
      }
   }

   printClassData(tbox::plog);

}

}
}
#endif // included_EmbeddedBoundaryShapePolygon_C
