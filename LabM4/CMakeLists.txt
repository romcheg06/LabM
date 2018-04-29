set(APP_NAME labM4)

set(HEADERS ../Utils/utils.h
            ../Utils/slice.h
            ../Utils/extendedslice.h
            lab4types.h)

FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
link_directories(${Boost_LIBRARY_DIR})

add_executable(${APP_NAME} main.cpp ${HEADERS})
target_link_libraries(${APP_NAME} ${MPI_CXX_LIBRARIES} stdc++ ${Boost_LIBRARIES})

if(MPI_CXX_COMPILE_FLAGS)
  set_target_properties(${APP_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
  set_target_properties(${APP_NAME} PROPERTIES
    LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()


