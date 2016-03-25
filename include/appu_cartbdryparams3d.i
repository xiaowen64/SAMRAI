c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/apputils/boundary/appu_cartbdryparams3d.i $
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for 3d boundary constant common blocks
c
      common/cartbdrylocparams3d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT, 
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1  
c
c
      common/cartbdrycondparams3d/
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
