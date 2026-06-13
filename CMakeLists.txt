cmake_minimum_required(VERSION 3.16)

project(laba_pro LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(laba_pro
    main.cpp
    DataStructs.h
    GuidanceMethods.h
    GuidanceMethods.cpp
    KinematicFunc2D.h
    KinematicFunc2D.cpp
    LibConstFunc.h
    CalcTrajectory.h
    CalcTrajectory.cpp
)

include(GNUInstallDirs)
install(TARGETS laba_pro
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
