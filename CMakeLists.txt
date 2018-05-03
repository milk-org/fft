project(lib_fft_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")

# PROJECT_SOURCE_DIR
add_library(fft SHARED fft.c)
set_target_properties(fft PROPERTIES COMPILE_FLAGS "-DFFTCONFIGDIR=\\\"${PROJECT_SOURCE_DIR}/config\\\"")

install(TARGETS fft DESTINATION lib)
install(FILES fft.h DESTINATION include)

