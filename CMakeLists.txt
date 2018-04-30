project(lib_fft_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")
add_library(fft SHARED fft.c)

install(TARGETS fft DESTINATION lib)
install(FILES fft.h DESTINATION include)

