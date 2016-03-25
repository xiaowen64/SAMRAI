//
//
// This class creates mapped multi-block grid geometrys. The supported grid types
// include Cartesian, Wedge, and Spherical shell.  The spherical shell
// case is a full multi-block grid with 3 blocks.
//
// 

#ifndef included_MblkGeometryXD
#define included_MblkGeometryXD

#include "SAMRAI_config.h"

#include "tbox/Array.h"
#include "Box.h"
#include "BoxArray.h"
#include "CellData.h"
#include "tbox/Database.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"
#include "BlockGridGeometry.h"

using namespace SAMRAI;

class MblkGeometry
{
public:

   //
   // Reads geometry information from the "MblkGeometry" input file
   // entry.                                                                 
   //
   MblkGeometry(const string& object_name,
                tbox::Pointer<tbox::Database> input_db,
                const int nblocks);
   
   ~MblkGeometry();
   
   //
   // Return the geometry type (CARTESIAN, WEDGE, or SPHERICAL_SHELL)
   //
   string getGeometryType();

   //
   // Return the user-specified refine boxes, given a block and 
   // level number
   //
   bool getRefineBoxes(hier::BoxArray<NDIM>& refine_boxes,
                       const int block_number,
                       const int level_number);
   
   //
   // Build mapped grid on patch.  The method defers the actual grid 
   // construction to private members, depending on the geometry 
   // choice in input.
   //
   void buildGridOnPatch(const hier::Patch<NDIM>& patch,
                         const hier::Box<NDIM>& domain,
                         const int xyz_id,
                         const int level_number,
                         const int block_number,
			 int *dom_local_blocks );

   //
   // The array of block indices denoting patch neighbors
   //
   void buildLocalBlocks( const hier::Box<NDIM>& pbox,   // the patch box
			  const hier::Box<NDIM>& domain, // the block box
			  const int block_number,
			  int *dom_local_blocks );  // this returns the blocks neighboring this patch

private:

   //
   // Read data members from input.  
   //
   void getFromInput(tbox::Pointer<tbox::Database> input_db,
                     bool is_from_restart);
   
   //
   // the cartesian input
   //
   void buildCartesianGridOnPatch(const hier::Patch<NDIM>& patch,
				  const hier::Box<NDIM>& domain,
                                  const int xyz_id,
                                  const int block_number);

   
   //
   // Wedge grid construction.
   //
   void buildWedgeGridOnPatch(const hier::Patch<NDIM>& patch,
			      const hier::Box<NDIM>& domain,
                              const int xyz_id,
                              const int block_number);


   //
   // trilinearly interpolated base mesh
   //
   void buildTrilinearGridOnPatch(const hier::Patch<NDIM>& patch,
				  const hier::Box<NDIM>& domain,
				  const int xyz_id,
				  const int block_number );

   //
   // Spherical shell grid construction
   //
   void buildSShellGridOnPatch(const hier::Patch<NDIM>& patch,
                               const hier::Box<NDIM>& domain,
                               const int xyz_id,
			       const int block_number,
			       const int level_number );

   //
   // For the spherical shell construction, i always points in the r direction
   // and j,k are points on the shell.  Given a certain j,k compute the 
   // unit sphere locations for block face (actual xyz is computed
   // by x = r*xface, y = r*yface, z = r*zface.  Note that the dimension
   // in the theta direction (nth) should be the same for each block.
   //
   void computeUnitSphereOctant(int nblock, 
                                int nth, 
                                int j, int k, 
                                double *xface, 
                                double *yface, 
                                double *zface );

   //
   // Geometry type.  Choices are CARTESIAN, WEDGE, SPHERICAL_SHELL
   //
   string d_geom_problem;
   string d_object_name;  // name of the object to pull in data from input parser
   

   //
   // The number of blocks and the set of skelton grid geometries that make
   // up a multiblock mesh.
   //
   int d_nblocks;
   
   //
   // Cartesian inputs
   //
   tbox::Array<tbox::Array<double> > d_cart_xlo;
   tbox::Array<tbox::Array<double> > d_cart_xhi;

   //
   // Wedge inputs
   //
   tbox::Array<double> d_wedge_rmin;
   tbox::Array<double> d_wedge_rmax;
   double d_wedge_thmin;
   double d_wedge_thmax;
   double d_wedge_zmin;
   double d_wedge_zmax;

   //
   // trilinear inputs
   //
   string d_tri_mesh_filename;
   int d_tri_nblocks;      // number of blocks
   int *d_tri_nxp;          // block bounds
   int *d_tri_nyp;
   int *d_tri_nzp;
   int *d_tri_node_size;   // block size

   int **d_tri_nbr;   // integer array of neighboring blocks
   double **d_tri_x; // [block][node]  nodal coordinates
   double **d_tri_y;
   double **d_tri_z;
   
   //
   // Shell inputs
   //
   double d_sshell_rmin;
   double d_sshell_rmax;

   // options are SOLID, OCTANT
   string d_sshell_type;

   //
   // For tagging in the spherical octant case
   //
   double d_tag_velocity;
   double d_tag_width;


   // if SOLID, need to read in these...
   double d_sangle_degrees;
   double d_sangle_thmin;
   
   //
   // Refine boxes for different blocks/levels
   //
   tbox::Array< tbox::Array< hier::BoxArray<NDIM> > > d_refine_boxes;
   
};


#endif // included_MblkGeometry
