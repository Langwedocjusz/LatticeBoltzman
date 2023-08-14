cmake_minimum_required(VERSION 3.16.3)

#Set the build output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#Print selected build type / available build types in case of Visual Studio
if (CMAKE_GENERATOR MATCHES "Visual Studio")
	message(STATUS "Build types: ${CMAKE_CONFIGURATION_TYPES}")
else()
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

#Set the project name (and version)
project(LatticeBoltzman)

#Specify CPP standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Build current project as an executable
add_executable(${PROJECT_NAME})

#Specify source files
file(GLOB_RECURSE headers src/*.h)
file(GLOB_RECURSE sources src/*.cpp)

target_sources(${PROJECT_NAME} PRIVATE ${headers} ${sources} ${imgui_impl})

#Specify include directories
target_include_directories(${PROJECT_NAME} PUBLIC src)

#Add nlohmann json
target_include_directories(${PROJECT_NAME} PUBLIC vendor/json/single_include)

#Directory structure for IDEs like Visual Studio
source_group(src REGULAR_EXPRESSION "src/*")

#Set current project as Visual Studio start project
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
#Set Visual Studio Debugger working directory
set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

#Need to link agains tbb to make code using parallel for_each work
#with older g++ versions
if (NOT WIN32)
	find_package(TBB REQUIRED)
	target_link_libraries(${PROJECT_NAME} PRIVATE TBB::tbb)
endif()
