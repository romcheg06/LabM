set(LAB_M2_NAME labM2)

set(HEADERS ../Utils/utils.h slice.h extendedslice.h)

FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
link_directories(${Boost_LIBRARY_DIR})

add_executable(${LAB_M2_NAME} main.cpp ${HEADERS})
target_link_libraries(${LAB_M2_NAME} ${MPI_CXX_LIBRARIES} stdc++ ${Boost_LIBRARIES})

if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(${LAB_M2_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(${LAB_M2_NAME} PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()

add_subdirectory(FieldViewer)

