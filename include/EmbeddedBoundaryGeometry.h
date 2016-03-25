//
// File:        EmbeddedBoundaryGeometry.h
// Package:     SAMRAI 
//              Structured Adaptive Mesh Refinement Applications Infrastructure
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Release:     $Name:  $
// Revision:    $Revision: 609 $
// Modified:    $Date: 2005-09-13 15:15:49 -0700 (Tue, 13 Sep 2005) $
// Description: Construction and management of embedded boundary data 
//              on an AMR hierarchy.
// 

#ifndef included_appu_EmbeddedBoundaryGeometry
#define included_appu_EmbeddedBoundaryGeometry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_appu_CutCell
#include "CutCell.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_geom_CartesianGridGeometry
#include "CartesianGridGeometry.h"
#endif
#ifndef included_pdat_CellVariable
#include "CellVariable.h"
#endif
#ifndef included_appu_CubesPatchInterface
#include "CubesPatchInterface.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_appu_EmbeddedBoundaryDefines
#include "EmbeddedBoundaryDefines.h"
#endif
#ifndef included_appu_EmbeddedBoundaryShape
#include "EmbeddedBoundaryShape.h"
#endif
#ifndef included_appu_ElevenPatchInterface
#include "ElevenPatchInterface.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_pdat_IndexVariable
#include "IndexVariable.h"
#endif
#ifndef included_pdat_NodeData
#include "NodeData.h"
#endif
#ifndef included_pdat_NodeVariable
#include "NodeVariable.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_xfer_RefineAlgorithm
#include "RefineAlgorithm.h"
#endif
#ifndef included_xfer_RefinePatchStrategy
#include "RefinePatchStrategy.h"
#endif
#ifndef included_tbox_Serializable
#include "tbox/Serializable.h"
#endif
#ifndef included_tbox_Timer
#include "tbox/Timer.h"
#endif
#ifndef included_appu_VisItDataWriter
#include "VisItDataWriter.h"
#endif
#ifndef included_appu_VisMaterialsDataStrategy
#include "VisMaterialsDataStrategy.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
   namespace appu {

/*!
 * @brief Class EmbeddedBoundaryGeometry provides embedded boundary mesh
 * construction, storage, and management on an AMR hierarchy.  
 * 
 * The embedded boundary may be constructed from a set of analytic shapes 
 * supplied through input.  The following outlines the steps required:
 *
 * 1. Construct an EmbeddedBoundaryGeometry object, supplying the input file
 *    that contains the shape entries, a pointer to the Cartesian grid 
 *    geometry, and the desired number of ghosts to be used in defining
 *    the embedded boundary:
 *
 *    \verbatim
 *    EmbeddedBoundaryGeometry* eb_geom = 
 *       new EmbeddedBoundaryGeometry("EmbeddedBoundaryGeometry",
 *                                     input_db->getDatabase("EBdryGeometry"),
 *                                     grid_geometry,
 *                                     nghosts);
 *    \endverbatim
 *
 *    Note:  The cart_grid_geometry argument is optional and may be supplied 
 *           later using the "setGridGeometry()" method.  However, it must
 *           be set before the "buildEmbeddedBoundaryOnLevel()" is called.
 *           The nghosts argument is also optional and is set to zero by
 *           default.
 *
 * 2. Build the embedded boundary on the levels of the hierarchy:
 *
 *    \verbatim
 *    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
 *    eb_geom->buildEmbeddedBoundaryOnLevel(level);
 *    \endverbatim
 *
 *    Note: It is also possible to pass in a hierarchy and an old level
 *          as arguments.  If supplied, this information can be used to 
 *          accelerate construction of the boundary on the supplied level.
 *
 * 3. Access information about the embedded boundary from patches:
 *   
 *    - The "cell flag" identifies whether a cell is cut, solid, or flow
 *    - The "node flag" identifies nodes as inside, outside, boundary 
 *      (first node just inside solid boundary), or on-boundary (within
 *      a specified distance from the solid boundary).
 *    - The "volume fraction" will be 0.0 for solid cells, 1.0 for flow
 *      cells, and somewhere in-between for cut cells
 *    - The list of "cut cells" holds information about the normal, face 
 *      areas, etc. on specific cells that are cut.
 *
 *    \verbatim
 *    int cell_flag_index = eb_geom->getCellFlagDataId();  
 *    int node_flag_index = eb_geom->getNodeInsideOutsideDataId();  
 *    int vol_frac_index = eb_geom->getCellVolumeDataId();  
 *    int cut_cell_index = eb_geom->getIndexCutCellDataId();
 *
 *    tbox::Pointer<CellData<NDIM,int> > cell_flag_data = 
 *        patch->getPatchData(cell_flag_index); 
 *    tbox::Pointer<NodeData<NDIM,int> > node_flag_data = 
 *        patch->getPatchData(node_flag_index); 
 *    tbox::Pointer<CellData<NDIM,double> > vol_frac_data = 
 *        patch->getPatchData(vol_frac_index); 
 *    tbox::Pointer<IndexData<NDIM,CutCell> > cut_cell_data = 
 *        patch->getPatchData(cut_cell_index); 
 *    \endverbatim
 *
 * Input parameters specify the list of shapes to be used to construct the 
 * boundary, the desired accuracy of the volume and area fraction 
 * computation, and information about whether write constructed embedded
 * boundary information to file, or read from file.
 *
 * Required input keys and data types: NONE
 *
 * Optional input keys, data types, and defaults:
 *
 *    -\b verbose
 *      boolean specifying whether to output information about the
 *      embedded boundary, such as the number of cut cells, the
 *      error in the volume calculation, etc. to the
 *      log file. If no value is supplied in input, the default
 *      is TRUE.
 *
 *    -\b max_subdivides  
 *      integer specifying the number of cell subdivides when
 *      computing the volume and area fractions.  The larger
 *      the number of subdivides, the more accurate the fraction 
 *      calculation will be, but the calculation will also be more 
 *      expensive.  If no value is supplied in input, a default of
 *      0 is used.
 *
 *    -\b read_from_file 
 *      bool specifying whether to read the embedded boundary from file.  
 *      If true, it will read flag, volume fraction, and cut cell information
 *      describing the embedded boundary from the specified HDF "dirname" 
 *      (specified below). If no value is supplied, the default is FALSE. 
 *
 *    -\b write_to_file 
 *      bool specifying whether to write the embedded boundary to the
 *      specified HDF "dirname" (specified below).  If true, it will write 
 *      flag, volume fraction, and cut cell information computed while 
 *      building the embedded boundary to the specified HDF "dirname" 
 *      (specified below). If no value is supplied, the default is FALSE. 
 *
 *    -\b dirname
 *      string specifying the name of the HDF directory to read/write
 *      embedded boundary information.  Used with the "read_from_file" 
 *      and "write_to_file" options.  
 *
 *    -\b compute_areas_and_normal  
 *      bool specifying whether to compute areas and normal information.  
 *      Some applications only need the volume fraction so it is not 
 *      necessary to invoke the extra cost of computing the areas and 
 *      normal.  If no value is supplied, the default is TRUE.
 *
 *    -\b compute_cutcell_index_data 
 *      bool specifying whether to compute a list of cut cells.  If false, 
 *      it will compute the cell flag and volume fraction information but 
 *      will not create and store the list of CutCell data structs.  If
 *      no value is supplied, the default is TRUE. 
 *
 *    -\b compute_boundary_node_data
 *      bool specifying whether to mark nodes that are just inside the 
 *      geometry as BOUNDARY nodes, and to compute the centroid of the
 *      wetted cut area (i.e. the "front area") of the cut cells.  This 
 *      information may be needed for node-based finite difference or 
 *      finite element computations of the embedded boundary. If
 *      no value is supplied, the default is TRUE. 
 *
 *    -\b use_recursive_algs 
 *      bool specifying whether to use a recursive algorithm to compute 
 *      volume and area fractions.  If true, it volume and area fractions
 *      are computed by recursively subdividing the cell until the max 
 *      number of subdivides is reached.  If false, it will apply the max
 *      subdivides to divide the cell into a (potentially large) array of 
 *      subcells.  Algorithmically, the recursive algorithm has fewer 
 *      operations but the non-recursive algorithm may be more computationally 
 *      efficient because it can do array-based operations.  If
 *      no value is supplied, the default is TRUE. 
 *
 *    -\b Shapes 
 *      sub-database that specifies information about the analytic
 *      shapes used to construct an embedded boundary.  See the 
 *      EmbeddedBoundaryShapeSphere and EmbeddedBoundaryShapePolygon class
 *      headers for information on the inputs required.
 *
 *    -\b CubesPatchInterface
 *      sub-database for the Cubes interface, a cut-cell mesh generator 
 *      from NASA Ames.  See the CubesPatchInterface class header for 
 *      information about the required inputs.
 *
 *       
 * The following represents a sample input entry: 
 * 
 * \verbatim
 *  EmbeddedBoundaryGeometry{
 *     max_subdivides = 2
 *     read_from_file  = FALSE 
 *     write_to_file   = FALSE
 *     dirname         = "eb_grid"
 *     compute_areas_and_normal   = TRUE 
 *     compute_boundary_node_data = FALSE
 *     use_recursive_algs         = FALSE 
 *
 *     Shapes {
 *        Shape1 {
 *           type = "POLYGON"    
 *           vertices {
 *              v1 = 1.0 , 1.0
 *              v2 = 1.5 , .5 
 *              v3 = 2.0 , 1.75 
 *              v4 = 1.5 , 4.5
 *              v5 = .5 , 2.0 
 *           }
 *           height = 5.0 // only used for 3D 
 *        }
 *        Shape2 {
 *           type = "SPHERE"    
 *           center = 65., 50.
 *           radius = 20.     
 *        }
 *     }
 *  }
 *  \endverbatim
 *
 * Note: Each shape has its own specific set of inputs.  See the 
 *       individual shapes for information about input requirements.
 *
 * @see appu::EmbeddedBoundaryShapeSphere 
 * @see appu::EmbeddedBoundaryShapePolygon 
 * @see appu::CubesPatchInterface 
 * @see appu::ElevenPatchInterface 
 * @see appu::CutCell 
 */

template<int DIM> class EmbeddedBoundaryGeometry 
:
public appu::VisMaterialsDataStrategy<DIM>,
public xfer::RefinePatchStrategy<DIM>,
public tbox::Serializable
{
public:

   /*!
    * Enumerated type for the different cell classifications.
    *
    * - \b SOLID          {Cell is located in the "solid" region.}
    * - \b CUT            {Cell is cut, meaning a CutCell data
    *                      structure will be maintained at this cell.}
    * - \b BORDER         {Cell neighbors a cut cell, in the "flow" region.}
    * - \b FLOW           {Cell is located in the "flow" region.}
    */
   enum CELL_TYPE{ SOLID  = EmbeddedBoundaryDefines::SOLID,
                   CUT    = EmbeddedBoundaryDefines::CUT,
                   BORDER = EmbeddedBoundaryDefines::BORDER,
                   FLOW   = EmbeddedBoundaryDefines::FLOW };

  /*!
    * Enumerated type for inside/outside node classification.
    *
    * - \b INSIDE         {Node is located "inside" the prescribed geometry.}
    * - \b OUTSIDE        {Node is outside the geometry.}
    * - \b BOUNDARY       {Node is on the boundary of the geometry.  That is
    *                      it is the first one "inside" the geometry.}
    * - \b ONBOUNDARY     {Node is located exactly on the boundary of the
    *                      geometry (used to avoid divide-by-zero problems)
    *                      in numerical operations at embedded boundary.}
    */
   enum NODE_TYPE{ OUTSIDE    = EmbeddedBoundaryDefines::OUTSIDE,
                   INSIDE     = EmbeddedBoundaryDefines::INSIDE,
                   BOUNDARY   = EmbeddedBoundaryDefines::BOUNDARY,
                   ONBOUNDARY = EmbeddedBoundaryDefines::ONBOUNDARY};     

   /*!
    * Constructor sets default values and reads data from input.
    *
    * @param object_name  Name of object.
    * @param input_db     Input database.
    * @param grid_geom    The grid geometry (e.g. cartesian) used 
    *                     in the problem.
    * @param nghosts      Number of ghosts used to hold the embedded
    *                     boundary.  If not supplied, defaults to 0.
    *
    * The grid_geom and nghosts may be NULL if the embedded
    * boundary information is to be supplied by a file, either by restart
    * or by other input such as CART3D.  If the embedded boundary is to
    * be constructed using analytic shapes in SAMRAI, the input_db and
    * and grid_geom arguments must be supplied and may not be NULL.
    */
   EmbeddedBoundaryGeometry(
      const string& object_name,
      tbox::Pointer<tbox::Database> input_db = NULL,
      const tbox::Pointer<geom::CartesianGridGeometry<DIM> > grid_geom = NULL,
      const hier::IntVector<DIM>& nghosts = hier::IntVector<DIM>(0));


   /*!
    * Destructor deallocates data describing and unregisters the 
    * object with the restart manager if previously registered.
    */
   ~EmbeddedBoundaryGeometry();
   
   /*!
    * Build an embedded boundary by forming the set of cut cells, extend 
    * them to appropriately apply the physical boundary conditions, and 
    * lastly compute the surrounding volumes for mass correction.  
    *
    * Depending on the arguments supplied, this method may be used in 
    * one of two ways.  The first, taking only the level as an argument,
    * performs an exhaustive search of all cells on the level
    * to find the cut cells and classify the cells as inside or outside.  
    * The second, taking as additional arguments a hierarchy, coarser_level,
    * and possibly old_level, performs the same function but uses the 
    * information on the coarser and old levels to narrow the search for 
    * cut cells, making it considerably faster.  Generally, the first method 
    * is used for the coarsest level only and the second is used for all 
    * subsequent finer levels.             
    *
    * @param level Patch level on which embedded boundary is to be 
    *              constructed.
    * @param hierarchy Patch hierarchy of the level.  Required if the 
    *                  level supplied in the first argument is in the 
    *                  hierarchy and is not the coarsest level.
    * @param old_level Patch level which holds "old" embedded boundary
    *                  data (not required).  Use this if regridding
    *                  and you have an embedded boundary at the old
    *                  level to speed construction of the boundary on 
    *                  the new level.
    */
   void buildEmbeddedBoundaryOnLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      const tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy = NULL,
      const tbox::Pointer<hier::PatchLevel<DIM> > old_level = NULL);

   /*!
    * Tag nodes as being inside or outside the geometry.  Some applications
    * only require this knowledge and do not use need the volume/area
    * fraction information for cut cells. 
    *
    * @param level Patch level where tagging takes place. 
    */
   void tagInsideOutsideNodesOnLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > level);

   /*!
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with the VisIt 
    * visualization tool.
    *
    * @param visit_writer    VisIt data writer
    */
#ifdef HAVE_HDF5
   void registerVisItDataWriter( 
      tbox::Pointer<appu::VisItDataWriter<DIM> > visit_writer);
#endif

   /*!
    * Return the descriptor index for the CellData<DIM,int> patch data 
    * that holds the integer flag at each cell, identifying it as 
    * being solid, cut, boundary, or flow.  
    */
   int getCellFlagDataId() const;

   /*!
    * Return the descriptor index for the CellData<DIM,double> patch data 
    * that holds the cell volume fraction.
    */
   int getCellVolumeDataId() const;

   /*!
    * Return the descriptor index for the IndexData<DIM,CutCell> patch
    * data that holds the list of embedded boundary cells on the patch.
    */
   int getIndexCutCellDataId() const;

   /*!
    * Return the descriptor index for the NodeData<DIM,int> patch
    * data that specifies a cell as being inside or outside the geometry.
    */
   int getNodeInsideOutsideDataId() const;

   /*!
    * Compute the total volume, which is the sum of the volume fractions,
    * on the supplied level.
    */
   double computeTotalVolumeOnLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > level);

   /*!
    * Set the grid geometry, if it was not supplied when the object
    * was constructed.
    */
   void setGridGeometry(
      const tbox::Pointer<geom::CartesianGridGeometry<DIM> > grid_geom);
   
   /*!
    * Write embedded boundary information - cell flag and cut cell
    * information - to supplied directory.   Files of the form 
    * "ebmesh-l<ln>-p<pid>.hdf", where ln is the level number
    * and pid is the processor number, will be written to the directory.
    */
   void writeLevelEmbeddedBoundaryDataToFile(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      const string& dirname) const;
 
   ///
   ///  The following routine:
   ///
   ///      packMaterialFractionsIntoDoubleBuffer()
   ///
   ///  is a concrete implementation of function declared in the
   ///  appu::VisMaterialsDataStrategy abstract base class.
   ///
   
   /*!
    * Put volume data located on the patch into the double buffer over
    * the specified region.
    *
    * @param dbuffer   double buffer into which materials data is packed.
    * @param patch     supplied patch on which materials data is defined
    * @param region    region over which data is packed
    * @param material_name    name of the material
    */
   int packMaterialFractionsIntoDoubleBuffer(
      double* dbuffer,
      const hier::Patch<DIM>& patch,
      const hier::Box<DIM>& region,
      const string& material_name);

   ///
   ///  The following routines:
   ///
   ///      setPhysicalBoundaryConditions(),
   ///      getRefineOpStencilWidth(),
   ///      postprocessRefine()
   ///
   ///  are concrete implementations of functions declared in the
   ///  xfer::RefinePatchStrategy abstract base class.
   ///
   
   /*!
    * Set the data in ghost cells corresponding to physical boundary
    * conditions.  Specific boundary conditions are determined by
    * information specified in input file and numerical routines.
    */
   void setPhysicalBoundaryConditions(hier::Patch<DIM>& patch,
                                      const double fill_time,
                                      const hier::IntVector<DIM>&
                                      ghost_width_to_fill);

   /*!
    * Perform user-defined refining operations.  This member function
    * is called before the other refining operators.  For this class,
    * no preprocessing is needed for the refine operators so it is 
    * setup to do nothing.
    */
   virtual void preprocessRefine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio)
      {
         (void) fine;
         (void) coarse;
         (void) fine_box;
         (void) ratio;
      }

   /*!
    * Postprocess data after the refinement operator is applied.  For
    * the embedded boundary data, refining involves two steps;  First,
    * refine the "flag" data, which will designate on the fine level
    * where the boundary exists;  Second, on each fine cell that is flagged
    * to possibly contain the boundary, go through and determine whether 
    * the cell is indeed cut by the boundary and adjust the volume, flag,
    * and boundary cell data on the fine level as necessary.  This method
    * invokes step 2 of this refine operation. 
    */
   virtual void postprocessRefine(
      hier::Patch<DIM>& fine,
      const hier::Patch<DIM>& coarse,
      const hier::Box<DIM>& fine_box,
      const hier::IntVector<DIM>& ratio);

   /*!
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual hier::IntVector<DIM> getRefineOpStencilWidth() const;

   ///
   ///  The following routine:
   ///
   ///      putToDatabase()
   ///
   ///  is a concrete implementation of a function declared in the
   ///  tbox::Serializable abstract base class.
   ///
   void putToDatabase(tbox::Pointer<tbox::Database> db);
   
private:

   /*
    * Initialize variables & communication schedules used to construct
    * the embedded boundary.  The supplied number of ghosts specifies
    * the ghost border over which the embedded boundary will be stored.  
    * This can vary depending on the needs of the application.  The
    * variables are registered with the supplied Visit or Vizamrai
    * data writer(s).
    *
    * @param nghosts Size of the ghost layer over which the embedded
    *                boundary will be stored. 
    */
   void initializeVariables(const hier::IntVector<DIM>& nghosts);

   /*
    * Use CUBES or ELEVEN packages to construct the embedded boundary 
    * on the level.
    */
   void computeEmbeddedBoundaryOnLevelWithPackage(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      int &cut_cells_on_level);

   /*
    * Use internal SAMRAI routines to generate the embedded boundary 
    * on the supplied level using an exhaustive search to classify 
    * every cell on the level.  
    */
   void computeEmbeddedBoundaryOnLevel(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      int &cut_cells_on_level,
      double &l2_volume_error,
      double &l2_area_error,
      double &max_volume_error,
      double &max_area_error);

   /*
    * Compute volume and area fraction, normals, etc. on the supplied
    * cut cell.
    */
   void calculateCutCellInformation(
      appu::CutCell<DIM>& cut_cell,
      const double *lower,
      const double *upper,
      const double& fullcellvol,
      const double *fullareas,
      double& volume_error_estimate,
      double& area_error_estimate);

   /*
    * Determine boundary nodes (those that lie immediately inside the
    * boundary) and compute the centroid of the wetted surface of the 
    * cut (i.e. the "front area") of the cut cell.
    */
   void calculateBoundaryNodeInformation(
      appu::CutCell<DIM>& cut_cell,
      tbox::Pointer< pdat::NodeData<DIM,int> >& node_flag,
      const double *lower,
      const double *upper,
      const double *dx);
   
   /*
    * Calculate the volume of cut-cell using a non-recursive algorithm.  
    * This method uses the same algorithm as the recursiveCalculateVolume()
    * method but is implemented in a non-recursive way.  In practice, it
    * uses slightly more memory but is computationally faster than the
    * recursive version.  
    *
    * The computed volume is returned and the error_estimate is a 
    * call-by-reference argument that is set inside the method.
    */
   double calculateVolume(const double *cell_lower,
                          const double *cell_upper,
                          double& error_estimate) const;  

   /*
    * Calculate the volume of the cut-cell using a recursive algorithm.  
    * The cell is recursively sub-divided into smaller cells to identify
    * the cut plane.  The volume and error estimate args are set inside 
    * the method.  The return value indicates if it has reached the 
    * maximum number of subdivides for the cell - return value is 0 if 
    * it has NOT reached max subdivides, 1 if it has.
    */
   int recursiveCalculateVolume(const double *cell_lower,
                                const double *cell_upper,
                                double& volume,
                                double& error_estimate,
                                int& subdivide_level) const;   

   /*
    * Calculate face area using a non-recursive algorithm. This method 
    * uses the same algorithm as the recursiveCalculateArea()
    * method but is implemented in a non-recursive way.  In practice, it
    * uses slightly more memory but is computationally faster than the
    * recursive version.  
    *
    * The computed area is returned and the error_estimate is a 
    * call-by-reference argument that is set inside the method.
    */
   double calculateArea(const double *cell_lower,
                        const double *cell_upper,
                        const int face,
                        double& error_estimate) const;

   /*
    * Calculate face area using a recursive algorithm.  The face is 
    * recursively subdivided into smaller faces  The area and error 
    * estimate args are set inside the method.  The return value indicates 
    * if it has reached the maximum number of subdivides for the cell.
    * The return value is 0 if it has NOT reached max subdivides, 1 if it has.
    */
   int recursiveCalculateArea(const double *cell_lower,
                              const double *cell_upper,
                              const int face,
                              double& area,
                              double& error_estimate,
                              int& subdivide_level) const;


   /*
    * Classify cell as being FLOW, SOLID, or CUT.
    */
   int classifyCell(const double *lower, 
                    const double *upper) const;
   
   /*
    * Invoke refine algorithm to refine the embedded boundary definition
    * from coarser.
    */
   void refineEmbeddedBoundary(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      const tbox::Pointer<hier::PatchLevel<DIM> > old_level = NULL,
      const tbox::Pointer<hier::PatchHierarchy<DIM> > hierarchy = NULL);

   /*!
    * Set the embedded boundary cells on physical boundarys according
    * to the specified boundary conditions.
    */
   void setEmbeddedBoundaryAtPhysicalBoundaries(
      const tbox::Pointer<hier::PatchLevel<DIM> > level);

   /*
    * Set the level ratio information needed to construct the embedded
    * boundary.  This method should be invoked each time a new level
    * is added to the hierarchy which may have an embedded boundary
    * on it.
    */
   void setLevelRatioInformation(
      const tbox::Pointer<hier::PatchLevel<DIM> > level);
   
   /*
    * Set node inside/outside array from SAMRAI analytic shapes.
    */
   void doNativeShapeInsideOutside(
      const int* nx,
      const double* dx,
      const double* origin,
      int* node_flag) const;

   /*!
    * Read embedded boundary information - cell flag and cut cell
    * information - from supplied directory.  Files of the form 
    * "ebmesh-l<ln>-p<pid>.hdf", where ln is the level number
    * and pid is the processor number, should exist in the directory.
    * If they don't exist, an error will be thrown.
    */
   void readLevelEmbeddedBoundaryDataFromFile(
      const tbox::Pointer<hier::PatchLevel<DIM> > level,
      const string& dirname) const;

   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    */
   void getFromInput(tbox::Pointer<tbox::Database> db,
                     const bool is_from_restart);

   void getFromRestart();

   /*
    * Name used for this object.
    */
   string d_object_name;

   /*
    * Grid geometry which embedded boundary is applied.
    */
   tbox::Pointer<geom::CartesianGridGeometry<DIM> > d_grid_geometry;
   
   /*
    * Visit and Vizamrai data writers.
    */
#ifdef HAVE_HDF5
   tbox::Pointer<appu::VisItDataWriter<DIM> > d_visit_writer;   
#endif

   /*
    * Whether to write information about the embedded boundary construction
    * to log file.
    */
   bool d_verbose; 

   /*
    * INTERNAL or EXTERNAL flow.
    */
   string d_flow_type;          
   
   /*
    * Number of ghosts used for the embedded boundary.
    */
   hier::IntVector<DIM> d_ebdry_nghosts;
 
   /*
    * Data members that defined the embedded Boundary.
    */
   tbox::Pointer< pdat::IndexVariable<DIM,appu::CutCell<DIM> > > d_ebdry_var;
   tbox::Pointer< pdat::CellVariable<DIM,int> > d_cell_flag_var;
   tbox::Pointer< pdat::CellVariable<DIM,double> > d_cell_vol_var;
   tbox::Pointer< pdat::NodeVariable<DIM,int> > d_node_flag_var;
   
   int d_ebdry_data_id;
   int d_ebdry_scratch_id;
   int d_cell_flag_data_id;
   int d_cell_flag_scratch_id;
   int d_cell_vol_data_id;
   int d_cell_vol_scratch_id;
   int d_node_flag_data_id;

   /*
    * Refine algorithm and schedule to manage communication of embedded
    * boundary data from coarse to fine levels.
    */
   tbox::Pointer<xfer::RefineAlgorithm<DIM> > d_ebdry_refine_alg;

   /*
    * These should be provided by the gridding alg at some point.  We store
    * them locally for the time being.
    */
   int d_max_levels;
   tbox::Array<hier::IntVector<DIM> > d_ratio_to_coarser;
   tbox::Array<hier::IntVector<DIM> > d_level_ratios;
   
   /*
    * Interfaces to outside packages that are used to define the
    * embedded boundary.
    */
   bool d_use_cubes;
   tbox::Pointer<appu::CubesPatchInterface<DIM> > d_cubes_interface;
   
   bool d_use_eleven_inside_outside;
   bool d_use_eleven_boundary_node;
   tbox::Pointer<ElevenPatchInterface<DIM> > d_eleven_interface;

   /*
    * Array of native SAMRAI shapes used to define embedded boundary.
    */
   tbox::Array< tbox::Pointer< appu::EmbeddedBoundaryShape<DIM> > > d_shapes;

   /*
    * Max number of subdivisions allowed for computing volume and area
    * fractions.
    */
   int d_max_subdivides;

   /*
    * Bools specifying if we read/write EB data to/from file, and 
    * the directory name where the files are stored.
    */
   bool d_read_ebdry;
   bool d_write_ebdry;
   string d_ebdry_dirname;

   /*
    * Bools specifying if we should use recursive or non-recursive algorithms
    * to compute volume and areas of the cut cells, and whether or not
    * to compute the face areas and normal of the cut cell (some applications
    * don't need this information.
    */
   bool d_use_recursive_algs;
   bool d_compute_areas_and_normal;
   bool d_compute_cutcell_index_data;
   bool d_compute_boundary_node_data;
   
   /*
    * Boundary condition information.
    */
   tbox::Array<int> d_node_bdry_cond;
   tbox::Array<int> d_edge_bdry_cond;
   tbox::Array<int> d_face_bdry_cond;

   /*
    * Timers interspersed throughout the class.
    */
   tbox::Pointer<tbox::Timer> t_compute_eb;   
   tbox::Pointer<tbox::Timer> t_calc_node_inout;   
   tbox::Pointer<tbox::Timer> t_calc_volume;
   tbox::Pointer<tbox::Timer> t_calc_area;
   tbox::Pointer<tbox::Timer> t_calc_boundary_node;
   tbox::Pointer<tbox::Timer> t_eb_phys_bdry;
   tbox::Pointer<tbox::Timer> t_read_write_eb;

   /*
    * Used to record step count, for timing and statistics reporting.
    */
   int d_step_count;
   
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EmbeddedBoundaryGeometry.C"
#endif

