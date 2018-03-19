cmake_minimum_required(VERSION 3.8)
project(TPD_signal_decomposition)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp storage_adaptors.h
        unit_tests.h lin_alg.h lin_alg.h factor_operations.h gaussian_operations.h
        interpolation.h quadrature.h tpd_sim.h probabilistic_graph.h optimization.h solvers.h)
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
