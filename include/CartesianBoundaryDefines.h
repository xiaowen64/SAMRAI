//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/boundary/CartesianBoundaryDefines.h $
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Defines for boundary condition integer constants
//

#ifndef included_appu_CartesianBoundaryDefines
#define included_appu_CartesianBoundaryDefines

/*
 * Definitions for boundary types in 1d, 2d, and 3d:
 */

//@{
//! @name Definitions for boundary types in 1d, 2d, and 3d:
#define FACE3D_BDRY_TYPE     (1)
#define EDGE3D_BDRY_TYPE     (2)
#define NODE3D_BDRY_TYPE     (3)

#define EDGE2D_BDRY_TYPE     (1)
#define NODE2D_BDRY_TYPE     (2)

#define NODE1D_BDRY_TYPE     (1)
//@}


/*
 * Definitions for boundary array sizes in 1d, 2d, or 3d:
 */

//@{
//! @name Definitions for boundary array sizes in 1d, 2d, or 3d:
#define NUM_1D_NODES  (2)

#define NUM_2D_EDGES  (4)
#define NUM_2D_NODES  (4)

#define NUM_3D_FACES  (6)
#define NUM_3D_EDGES  (12)
#define NUM_3D_NODES  (8)
//@}


/*
 * Definitions for Face, Edge, and Node boundary locations:
 *
 * Note that these definitions are used only for:
 * - Node boundary locations in 1d (XLO, XHI only), or
 * - Edge boundary locations in 2d (XLO, XHI, YLO, YHI only), or
 * - Face boundary locations in 3d.
 */

//@{
//! @name Definitions for Face, Edge, and Node boundary locations (see source code for more information):
#define XLO        (0)
#define XHI        (1)
#define YLO        (2)
#define YHI        (3)
#define ZLO        (4)
#define ZHI        (5)
//@}


/*
 * Definitions for Node boundary locations in 2d:
 */

//@{
//! @name Definitions for Node boundary locations in 2d:
#define XLO_YLO_2D        (0)
#define XHI_YLO_2D        (1)
#define XLO_YHI_2D        (2)
#define XHI_YHI_2D        (3)
//@}


/*
 * Definitions for Edge boundary locations in 3d:
 */

//@{
//! @name Definitions for Edge boundary locations in 3d: 
#define YLO_ZLO_3D        (0)
#define YHI_ZLO_3D        (1)
#define YLO_ZHI_3D        (2)
#define YHI_ZHI_3D        (3)
#define XLO_ZLO_3D        (4)
#define XLO_ZHI_3D        (5)
#define XHI_ZLO_3D        (6)
#define XHI_ZHI_3D        (7)
#define XLO_YLO_3D        (8)
#define XHI_YLO_3D        (9)
#define XLO_YHI_3D        (10)
#define XHI_YHI_3D        (11)
//@}


/*
 * Definitions for Node boundary locations in 3d:
 */

//@{
//! @name Definitions for Node boundary locations in 3d: 
#define XLO_YLO_ZLO        (0)
#define XHI_YLO_ZLO        (1)
#define XLO_YHI_ZLO        (2)
#define XHI_YHI_ZLO        (3)
#define XLO_YLO_ZHI        (4)
#define XHI_YLO_ZHI        (5)
#define XLO_YHI_ZHI        (6)
#define XHI_YHI_ZHI        (7)
//@}


/*
 * Definitions for Face, Edge, and Node boundary conditions:
 *
 * Note that these definitions are used only for:
 * - Node boundary conditions in 1d, or
 * - Edge boundary conditions in 2d, or
 * - Face boundary conditions in 3d.
 */

//@{
//! @name Definitions for Face, Edge, and Node boundary conditions (see source code for more information):
#define FLOW_BC            (90) 
#define REFLECT_BC         (91)
#define DIRICHLET_BC       (92)
#define NEUMANN_BC         (93)
//@}


/*
 * Definitions for Edge and Node boundary conditions:
 * 
 * Note that the following definitions are used only for:
 * - Node boundary conditions in 2d (X and Y cases only), or
 * - Edge and Node boundary conditions in 3d.
 */

//@{
//! @name Definitions for Edge and Node boundary conditions (see source code for more information):
#define XFLOW_BC           (900)
#define YFLOW_BC           (901)
#define ZFLOW_BC           (902)
#define XREFLECT_BC        (910)
#define YREFLECT_BC        (911)
#define ZREFLECT_BC        (912)
#define XDIRICHLET_BC      (920)
#define YDIRICHLET_BC      (921)
#define ZDIRICHLET_BC      (922)
#define XNEUMANN_BC        (930)
#define YNEUMANN_BC        (931)
#define ZNEUMANN_BC        (932)
//@}

#endif
