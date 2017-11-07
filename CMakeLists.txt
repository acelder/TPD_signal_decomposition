cmake_minimum_required(VERSION 3.8)
project(TPD_signal_decomposition)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp finite_difference_simulations.h microporous_plate_system.cpp microporous_sphere_system.cpp mesoporous_film_system.cpp surface_system.cpp simulate_experiment.cpp)
add_executable(TPD_signal_decomposition ${SOURCE_FILES})

# find boost
find_package(Boost 1.63 EXACT)

if(Boost_FOUND)
    message("Boost found.")
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(TPD_signal_decomposition ${Boost_LIBRARIES})
else()
    message("Boost not found.")
endif()
