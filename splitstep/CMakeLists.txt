cmake_minimum_required (VERSION 3.0)
project (SplitStep)

#SET(CMAKE_CXX_COMPILER g++)

#find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    message("Found OpenMP")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories("../matrix")
add_subdirectory("../matrix" "matrix")

include_directories("../settings_manager")
add_subdirectory("../settings_manager" "settings_manager")

add_subdirectory("wavefunction")

add_executable(SplitStep main.cc arithmetics.cc ${SOURCE})
target_link_libraries(SplitStep fftw3)
target_link_libraries(SplitStep fftw3_omp)
target_link_libraries(SplitStep Matrix)
target_link_libraries(SplitStep SettingsManager)
target_link_libraries(SplitStep Wavefunction)
