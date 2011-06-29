/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Polygon embedded boundary shape 
 *
 ************************************************************************/

namespace SAMRAI {
namespace appu {

SAMRAI_INLINE_KEYWORD
void EmbeddedBoundaryShapePolygon::crossProduct(
   double a[3],
   const double b[3],
   const double c[3]) const
{
   a[0] = b[1] * c[2] - c[1] * b[2];
   a[1] = b[2] * c[0] - c[2] * b[0];
   a[2] = b[0] * c[1] - c[0] * b[1];
}

SAMRAI_INLINE_KEYWORD
double EmbeddedBoundaryShapePolygon::dotProduct(
   const double a[3],
   const double b[3]) const
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

SAMRAI_INLINE_KEYWORD
bool EmbeddedBoundaryShapePolygon::sameSide(
   double p1[3],
   double p2[3],
   double a[3],
   double b[3]) const
{
   double cp1[3], cp2[3], d;
   double b_minus_a[3];
   double p1_minus_a[3];
   double p2_minus_a[3];
   double dot1, dot2;

   b_minus_a[0] = b[0] - a[0];
   b_minus_a[1] = b[1] - a[1];
   b_minus_a[2] = b[2] - a[2];

   p1_minus_a[0] = p1[0] - a[0];
   p1_minus_a[1] = p1[1] - a[1];
   p1_minus_a[2] = p1[2] - a[2];

   p2_minus_a[0] = p2[0] - a[0];
   p2_minus_a[1] = p2[1] - a[1];
   p2_minus_a[2] = p2[2] - a[2];

   this->crossProduct(cp1, b_minus_a, p1_minus_a);
   this->crossProduct(cp2, b_minus_a, p2_minus_a);
   d = this->dotProduct(cp1, cp2);
   dot1 = this->dotProduct(cp1, cp1);
   dot2 = this->dotProduct(cp2, cp2);

   if (d >= -1000.e0 * d_eps * (dot1 + dot2))
      return true;
   else
      return false;

}

SAMRAI_INLINE_KEYWORD
bool EmbeddedBoundaryShapePolygon::pointInPolygon(
   tbox::Array<double> v_x,
   tbox::Array<double> v_y,
   double p_x,
   double p_y) const
{
   double p[3];
   double a[3];
   double b[3];
   double c[3];

   p[0] = p_x;
   p[1] = p_y;
   p[2] = 0.;
   for (int i = 0; i < d_num_vertices; i++) {
      a[0] = v_x[i];
      a[1] = v_y[i];
      a[2] = 0.;

      int ip1 = (i + 1) % d_num_vertices;
      int ip2 = (i + 2) % d_num_vertices;
      b[0] = v_x[ip1];
      b[1] = v_y[ip1];
      b[2] = 0.;
      c[0] = v_x[ip2];
      c[1] = v_y[ip2];
      c[2] = 0.;

      if (!sameSide(p, a, b, c))
         return false;
   }
   return true;
}

SAMRAI_INLINE_KEYWORD
bool EmbeddedBoundaryShapePolygon::isInside(
   const double* xyz) const
{

   bool below_height = true;
   bool inside_poly = false;

   /*
    * Check if height z is less than d_height.  If not, no further checks
    * are necessary.
    */
   if (d_dim == tbox::Dimension(3)) {
      below_height = xyz[2] < d_height;
   }

   /*
    * See if supplied point is inside or outside the polygon.
    */
   if (below_height) {
      inside_poly = pointInPolygon(d_vx, d_vy, xyz[0], xyz[1]);
   }

   return inside_poly;

}

SAMRAI_INLINE_KEYWORD
void EmbeddedBoundaryShapePolygon::isInside(
   const int* nx,
   const double* dx,
   const double* origin,
   int* inout) const
{

   /*
    * See if a set of points are inside or outside the polygon.
    */
   int kstr = 0;
   int jstr = 0;
   int istr = 0;
   int ijk;
   double xyz[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   bool below_height;
   bool inside_poly;

   /*
    * Initialize xyz to avoid uninitialized use.
    */
   for(int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      xyz[d] = 0.0;
   }

   /*
    * Compute the node centered inout points
    */
   int khi = 1;
   if (d_dim == tbox::Dimension(3)) {
      khi = nx[2];
   }
   for (int kk = 0; kk < khi; kk++) {
      kstr = kk * nx[1] * nx[0];
      if (d_dim == tbox::Dimension(3)) {
         xyz[d_dim.getValue()
             - 1] = origin[d_dim.getValue() - 1] + dx[d_dim.getValue() - 1] * (double)kk;
      }
      for (int jj = 0; jj < nx[1]; jj++) {
         jstr = jj * nx[0];
         xyz[1] = origin[1] + dx[1] * (double)jj;
         for (int ii = 0; ii < nx[0]; ii++) {
            istr = ii;
            xyz[0] = origin[0] + dx[0] * (double)ii;

            below_height = true;
            inside_poly = false;

            if (d_dim == tbox::Dimension(3)) {
               below_height = xyz[d_dim.getValue() - 1] < d_height;
            }

            if (below_height) {
               inside_poly = pointInPolygon(d_vx, d_vy, xyz[0], xyz[1]);
            }

            ijk = kstr + jstr + istr;
            if (inside_poly) {
               inout[ijk] = EmbeddedBoundaryDefines::INSIDE;
               inout[ijk] = 1;
            }
         }
      }
   }
}

SAMRAI_INLINE_KEYWORD
const std::string& EmbeddedBoundaryShapePolygon::getObjectName() const
{
   return d_object_name;
}

} // namespace appu
} // namespace SAMRAI
