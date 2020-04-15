# library name
set(LIBNAME "milkfft")
set(SRCNAME "fft")


set(SOURCEFILES
	${SRCNAME}.c
	fft_autocorrelation.c
	fft_structure_function.c)

set(INCLUDEFILES
	${SRCNAME}.h)


# DEFAULT SETTINGS 
# Do not change unless needed
# =====================================================================

project(lib_${LIBNAME}_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")


# Library can be compiled from multiple source files
# Convention: the main souce file is named <libname>.c
# 
add_library(${LIBNAME} SHARED ${SOURCEFILES})

set_target_properties(${LIBNAME} PROPERTIES COMPILE_FLAGS "-DFFTCONFIGDIR=\\\"${PROJECT_SOURCE_DIR}/config\\\"")

install(TARGETS ${LIBNAME} DESTINATION lib)
install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})

