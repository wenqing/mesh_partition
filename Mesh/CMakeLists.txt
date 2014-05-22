# Source files

set( MESH_SOURCES
      Grain.h
      Node.h
	  Node.cpp
      Elem.h
      Elem.cpp
      Mesh.h
      Mesh.cpp
)

# Create the library
add_library(Mesh STATIC ${MESH_SOURCES} )

target_link_libraries(
	Mesh
)

