SOURCE_FILES ( ONCRPC_SRCS
bcopy.c
nt.c
Xdr_floa.c
Xdr_mem.c
Xdr_stdi.c
Xdr.c
Xdr_arra.c
Xdr_rec.c
Xdr_refe.c
Xdr_strarr.c
)

ADD_LIBRARY(oncrpc ONCRPC_SRCS)

IF (WIN32)
  LINK_LIBRARIES( wsock32 )
ENDIF (WIN32)

INSTALL_TARGETS(/lib/vtk oncrpc)



