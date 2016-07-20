# Source files
ADD_EXECUTABLE(${PROJECT_NAME} partition.cpp)
target_link_libraries(${PROJECT_NAME} Mesh)	

###
## Compile mpmetis
include_directories(${METIS_PATH}/libmetis)
include_directories(${METIS_PATH}/programs)
set(METIS_SOURCES
   ${METIS_PATH}/programs/mpmetis.c
   ${METIS_PATH}/programs/cmdline_mpmetis.c
   ${METIS_PATH}/programs/io.c
   ${METIS_PATH}/programs/stat.c
   )
add_executable(mpmetis ${METIS_SOURCES})
target_link_libraries(mpmetis metis)



