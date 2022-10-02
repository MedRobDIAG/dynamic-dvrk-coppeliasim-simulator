cmake_minimum_required(VERSION 3.5)
project(dvrkDynModelLib VERSION 0.1 LANGUAGES C CXX)

# CMake definitions
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(
-DNDEBUG
-DDVRKDYNMODELLIB_EXPORTS
-D_WINDOWS
-D_USRDLL
-D_CRT_SECURE_NO_WARNINGS
-DUSE_MATH_DEFINES
-DDO_NOT_USE_SHARED_MEMORY
-DNON_MATLAB_PARSING
-DMAX_EXT_API_CONNECTIONS=255
-D_CONSOLE
)

include_directories(
    $ENV{EIGEN}
    ${CMAKE_CURRENT_LIST_DIR}/dvrDynModelLib/
)

set(DYNLIB_HEADER_FILES
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDefines.hpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDynamics.hpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkKinematics.hpp
)

set(DYNLIB_SRCS_FILES
    ${DYNLIB_HEADER_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/pch.h
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/framework.h
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDynamics.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkKinematics.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/pch.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDynModelLib.cpp
)

set(TEST_SRCS_FILES
    ${SRCS_FILES}
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDefines.hpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDynamics.hpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkKinematics.hpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/pch.h
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/framework.h
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDynamics.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkKinematics.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/pch.cpp
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib/dvrkDynModelLib.cpp
)


add_library( dvrkDynModelLib STATIC
             ${DYNLIB_SRCS_FILES} )


# CMake instructions to test using the static lib

#set(COPPELIA '$ENV{VREPx64_coppelia}')
file(TO_CMAKE_PATH "$ENV{VREPx64_coppelia}" COPPELIA)

include_directories(
    ${COPPELIA}/programming/remoteApi
    ${COPPELIA}/programming/include
    ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib
)

add_executable( dvrkDynModelLib_test
                ${DYNLIB_HEADER_FILES}
                ${COPPELIA}/programming/remoteApi/extApi.h
                ${COPPELIA}/programming/remoteApi/extApiPlatform.h
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/utils.hpp
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/Timer.hpp
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/extApi.c
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/extApiPlatform.c
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/Timer.cpp
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/utils.cpp
                ${CMAKE_CURRENT_LIST_DIR}/dvrkDynModelLib_test/main.cpp
) 

target_link_libraries( dvrkDynModelLib_test
                       dvrkDynModelLib )
