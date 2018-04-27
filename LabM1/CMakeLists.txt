set(APP_NAME labM1)

set(HEADERS ../Utils/utils.h)

add_executable(${APP_NAME} main.cpp ${HEADERS})
target_link_libraries(${APP_NAME} ${MPI_CXX_LIBRARIES} stdc++)

if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(${APP_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(${APP_NAME} PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()

