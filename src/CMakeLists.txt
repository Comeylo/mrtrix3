## Version: $Id$
##
######################################################################
##
### Commentary:
##
######################################################################
##
### Change Log:
##
######################################################################
##
### Code:

add_subdirectory(core)
add_subdirectory(dwi)
add_subdirectory(gui)
add_subdirectory(stats)

## #################################################################
## Essential target
## #################################################################

set(MRTRIX3_TARGETS)

set(MRTRIX3_TARGETS ${MRTRIX3_TARGETS} mrtrix3Core)
set(MRTRIX3_TARGETS ${MRTRIX3_TARGETS} mrtrix3DWI)
set(MRTRIX3_TARGETS ${MRTRIX3_TARGETS} mrtrix3Gui)
set(MRTRIX3_TARGETS ${MRTRIX3_TARGETS} mrtrix3Stats)

## #################################################################
## Target export
## #################################################################

export(TARGETS ${MRTRIX3_TARGETS} FILE "${CMAKE_BINARY_DIR}/mrtrix3Depends.cmake")

######################################################################
### CMakeLists.txt ends here
