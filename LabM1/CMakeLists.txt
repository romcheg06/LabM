set(LAB_M1_NAME labM1)

set(HEADERS ../Utils/utils.h)

add_executable(${LAB_M1_NAME} main.cpp ${HEADERS})
target_link_libraries(${LAB_M1_NAME} ${MPI_CXX_LIBRARIES} stdc++)

if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(${LAB_M1_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(${LAB_M1_NAME} PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()

