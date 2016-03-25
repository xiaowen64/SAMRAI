INSTALL_PROGRAMS(/bin vizamrai vizamrai.1slice vizamrai.3slice vizamrai.carpet vizamrai.combo vizamrai.function vizamrai.iso vizamrai.vector vizamrai.vol vizamrai.tcl vizamrai.1slice.tcl vizamrai.3slice.tcl vizamrai.carpet.tcl vizamrai.combo.tcl vizamrai.function.tcl vizamrai.iso.tcl vizamrai.vector.tcl vizamrai.vol.tcl)


CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.1slice)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.3slice)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.carpet)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.combo)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.function)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.iso)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.vector)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/runvizamrai ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.vol)



CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.1slice.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.1slice.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.3slice.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.3slice.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.carpet.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.carpet.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.combo.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.combo.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.function.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.function.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.iso.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.iso.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.vector.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.vector.tcl)
CONFIGURE_FILE(${VIZAMRAI_SOURCE_DIR}/scripts/vizamrai.vol.tcl ${VIZAMRAI_BINARY_DIR}/scripts/vizamrai.vol.tcl)

