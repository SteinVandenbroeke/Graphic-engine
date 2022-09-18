cmake_minimum_required(VERSION 3.5)

project(engine)

############################################################
# Set compiler flags
############################################################
set(OWN_GXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS} -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
	easy_image.cc
	engine.cc
	ini_configuration.cc
	l_parser.cc
        L2.cpp L2.h
		vector3d.cc vector3d.h D3.cpp D3.h L3.cpp L3.h ZBuffer.cpp ZBuffer.h Light.cpp Light.h)


set( Qt5Core_DIR " /home/stein/Qt5.12.12/5.12.10/gcc_64/lib/cmake/Qt5Core ")
set( Qt5Widgets_DIR " /home/stein/Qt5.12.12/5.12.10/gcc_64/lib/cmake/Qt5Widgets ")
set( Qt5Gui_DIR " /home/stein/Qt5.12.12/5.12.10/gcc_64/lib/cmake/Qt5Gui ")
set( Qt5Multimedia_DIR " /home/stein/Qt5.12.12/5.12.10/gcc_64/lib/cmake/Qt5Multimedia ")

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )