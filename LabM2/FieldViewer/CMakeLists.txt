set(APP_NAME FieldViewer)

set(HEADERS mainwindow.h)
set(SOURCES main.cpp mainwindow.cpp)
set(FORMS mainwindow.ui)

include_directories(..)

FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
link_directories(${Boost_LIBRARY_DIR})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt5Widgets)

add_executable(${APP_NAME} ${SOURCES} ${HEADERS} ${FORMS})
target_link_libraries(${APP_NAME} ${Boost_LIBRARIES} Qt5::Widgets)


