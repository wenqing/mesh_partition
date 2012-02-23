# Source files


ADD_EXECUTABLE(${PROJECT_NAME} partition.cpp)


target_link_libraries(${PROJECT_NAME} 
	Math
	Mesh
)
