cmake_minimum_required(VERSION 3.10.0)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# List of .cpp files to be compiled
SET(SRCS
  mesh/halfedge_mesh.cpp
  mesh/meshio.cpp
  mesh/polygon_soup_mesh.cpp
  mesh/ply_halfedge_mesh_data.cpp
  mesh/mesh_embedded_curve.cpp
  
  geometry/geometry.cpp
  geometry/intrinsic_geometry.cpp
  geometry/geometry_cache.cpp 
  geometry/edge_length_geometry.cpp
  geometry/dependent_quantity.cpp

  numerical/suitesparse_utilities.cpp
  numerical/linear_solvers.cpp
  numerical/eigenproblem_solvers.cpp
  numerical/qr_solvers.cpp
  numerical/square_solvers.cpp
  numerical/positive_definite_solvers.cpp

  algorithms/mesh_graph_algorithms.cpp
  algorithms/discrete_operators.cpp
  algorithms/detect_symmetry.cpp
  algorithms/mesh_ray_tracer.cpp
  algorithms/direction_fields.cpp
  algorithms/exact_polyhedral_geodesics.cpp
  algorithms/fast_marching_method.cpp

  utilities/unit_vector3.cpp
  utilities/vector3.cpp
  utilities/vector2.cpp
  utilities/utilities.cpp
  utilities/quaternion.cpp
  utilities/disjoint_sets.cpp
)


# Create a single library for the project
add_library(
        geometry-central 
        ${SRCS}
        )

target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/eigen-git-mirror")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/nanort/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/nanoflann/include")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/glm")
target_include_directories(geometry-central PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/happly")
if(SUITESPARSE_FOUND)
  target_include_directories(geometry-central PUBLIC "${SUITESPARSE_INCLUDE_DIRS}")
endif()

# If we (optionally) resolved suitesparse support, we want to link against that as well
if(HAVE_SUITESPARSE)
    target_link_libraries(geometry-central ${SUITESPARSE_LIBRARIES})
endif(HAVE_SUITESPARSE)


# Set compiler properties for the library
set_property(TARGET geometry-central PROPERTY CXX_STANDARD 11)
