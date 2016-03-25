//
// File:        VisMaterialsDataStrategy.h
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2003 The Regents of the University of California
// Revision:    $Revision: 47 $
// Modified:    $Date: 2004-12-09 16:08:57 -0800 (Thu, 09 Dec 2004) $
// Description: Interface for writing material related data to a VisIt 
//              dump file.
//

#ifndef included_appu_VisMaterialsDataStrategy
#define included_appu_VisMaterialsDataStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_Patch
#include "Patch.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace appu {

/*!
 * @brief Class VisMaterialsDataStrategy<DIM> is an abstract base
 * class that defines an interface allowing an VisItDataWriter<DIM>
 * object to generate plot files that contain material and species
 * fractions, as well as state variables for individual materials.  A
 * concrete object of this type must be registered with the data
 * writer in order to use materials or species with the data writer.
 * The registration of the concrete object is done using the method
 * setMaterialsDataWriter() from the VisItDataWriter class.  VisIt
 * requires that material fractions, species fractions, and material
 * state variables be cell-centered.  If they are not cell-centered in
 * the simulation, it is the job of the relevant packing method to
 * convert them to a cell-centered basis before packing them into the
 * buffer.
 *
 * The concrete strategy object is responsible for supplying an
 * implementation of packMaterialFractionsIntoDoubleBuffer().  If
 * species are used, packSpeciesFractionsIntoDoubleBuffer() must also
 * be implemented.  If material state variables are used,
 * packMaterialStateVariableIntoDoubleBuffer() must be implemented.
 *
 * @see appu::VisItDataWriter
 */

template<int DIM> class VisMaterialsDataStrategy 
{
public:
   /*!
    * @brief Default constructor for VisMaterialsDataStrategy<DIM>.
    */
   VisMaterialsDataStrategy();

   /*!
    * @brief Destructor for VisMaterialsDataStrategy<DIM>.
    */
   virtual ~VisMaterialsDataStrategy<DIM>();

   /*!
    * @brief Enumerated type for the allowable return values
    * for packMaterialFractionsIntoDoubleBuffer() and
    * packSpeciesFractionsIntoDoubleBuffer().
    * - \b ALL_ZERO   (Fractions are 0.0 for every cell in patch.)
    * - \b ALL_ONE    (Fractions are 1.0 for every cell in patch.)
    * - \b MIXTURE    (There is some of this species/material in one or
    *                    more cells of this patch, but the above two cases
    *                   do not apply.)
    */

   enum PACK_RETURN_TYPE { VISIT_ALLZERO = 0,
                           VISIT_ALLONE = 1,
                           VISIT_MIXED = 2 };

   /*!
    * @brief This function, which must be implemented whenever
    * materials are used, packs cell-centered material fractions for
    * the given material, patch, and region into the given 1D double
    * precision buffer which will already have been allocated.  If a
    * non-zero ghost cell vector was specified when
    * registerMaterialNames() was invoked, then ghost data
    * corresponding to this ghost cell vector must be packed into this
    * double buffer.  The data must be packed into the buffer in
    * column major order, i.e. (f(x_0,y_0,z_0), f(x_1,y_0,z_0),
    * f(x_2,y_0,z_0), ...).
    *
    * This method will be called once for each material for each patch.
    *
    * A enumerated PACK_RETURN_TYPE is used for a return value.  To 
    * save space in the visit data file, you may choose to set the return
    * value to indicate if a material does not exist at all on the patch,
    * or if the material exists fully on the patch. A return of ALL_ZERO 
    * indicates there is 0% of the material in each of the cells of the 
    * patch, while a return type of ALL_ONE indicates the material consumes
    * 100% on each of the cells of the patch.  If the patch has a mixture
    * of the material (i.e. between 0% and 100%) then return MIXTURE.
    *
    * @param buffer Double precision array into which cell-centered
    *  material fractions are packed.
    * @param patch hier::Patch on which fractions are defined.
    * @param region hier::Box region over which to pack fractions.
    * @param material_name String identifier for the material.
    * @return The enumeration constant 
    *    VisMaterialsDataStrategy::ALL_ZERO, 
    *    VisMaterialsDataStrategy::ALL_ONE, 
    *    or VisMaterialsDataStrategy::MIXTURE.
    */
   virtual int packMaterialFractionsIntoDoubleBuffer(
      double *buffer,
      const hier::Patch<DIM>& patch,
      const hier::Box<DIM>& region,
      const string& material_name);

   /*!
    * @brief This function packs cell-centered species fractions for
    * the given species.
    *
    * This user supplied function packs species fractions of the 
    * given material, patch, and region into the supplied 1D double 
    * precision buffer.  If a non-zero ghost cell vector was specified when
    * registerSpeciesNames() was invoked, then ghost data
    * corresponding to this ghost cell vector must be packed into this
    * double buffer.  The data must be packed into the buffer in
    * column major order.
    *
    * This method will be called once for each species for each patch.
    *
    * The method must return a PACK_RETURN_TYPE of ALL_ONE, ALL_ZERO,
    * or MIXED.  See the discussion above for the 
    * "packMaterialFractionsIntoDoubleBuffer()" method for an explanation
    * of correct return values.
    *
    * @param buffer Double precision array into which  cell-centered 
    *   species fractions are packed.
    * @param patch hier::Patch on which fractions are defined.
    * @param region hier::Box region over which to pack fractions.
    * @param material_name String identifier for the material to 
    *  which the species belongs.
    * @param species_name String identifier for the species.
    * @return The enumeration constant 
    *    VisMaterialsDataStrategy::ALL_ZERO, 
    *    VisMaterialsDataStrategy::ALL_ONE, 
    *    or VisMaterialsDataStrategy::MIXTURE.
    */
   virtual int packSpeciesFractionsIntoDoubleBuffer(
      double *buffer,
      const hier::Patch<DIM>& patch,
      const hier::Box<DIM>& region,
      const string& material_name,
      const string& species_name);
};
}
}
#endif


#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "VisMaterialsDataStrategy.C"
#endif
