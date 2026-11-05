cmake_minimum_required(VERSION 3.28)
project(RDKitXTB VERSION 0.11.0 LANGUAGES CXX)
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(FATAL_ERROR "rdXTBTools currently supports Linux; the dynamic loader needs a platform implementation elsewhere")
endif()
include(CTest)
include(GNUInstallDirs)
find_package(RDKit CONFIG REQUIRED)
find_package(OpenSSL REQUIRED COMPONENTS Crypto)
find_path(XTB_INCLUDE_DIR xtb.h HINTS "${RDXTB_XTB_PREFIX}/include" REQUIRED)

add_library(RDKitXTB STATIC cpp/rdXTBTools/XTBTools.cpp)
set_target_properties(RDKitXTB PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(RDKitXTB PUBLIC cxx_std_17)
target_compile_options(RDKitXTB PRIVATE -Wall -Wextra)
target_include_directories(RDKitXTB PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/rdXTBTools>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rdkit-eth>)
target_include_directories(RDKitXTB PRIVATE "${XTB_INCLUDE_DIR}")
target_link_libraries(RDKitXTB PUBLIC RDKit::GraphMol OpenSSL::Crypto ${CMAKE_DL_LIBS})
if(RDXTB_XTB_PREFIX)
  target_compile_definitions(RDKitXTB PRIVATE RDXTB_DEFAULT_LIBRARY="${RDXTB_XTB_PREFIX}/lib/libxtb.so")
endif()

option(RDXTB_BUILD_PYTHON "Build the Boost.Python extension" ON)
if(RDXTB_BUILD_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
  find_package(rdkitpython CONFIG REQUIRED)
  Python_add_library(rdXTBTools MODULE WITH_SOABI cpp/rdXTBTools/Wrap.cpp)
  target_link_libraries(rdXTBTools PRIVATE RDKitXTB RDKit::RDBoost Python::NumPy)
  target_compile_definitions(rdXTBTools PRIVATE RDXTB_RDKIT_VERSION="${RDKit_VERSION}")
  target_compile_options(rdXTBTools PRIVATE -Wall -Wextra)
  set_target_properties(rdXTBTools PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

install(TARGETS RDKitXTB EXPORT RDKitXTBTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES cpp/rdXTBTools/XTBTools.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rdkit-eth)
install(EXPORT RDKitXTBTargets NAMESPACE RDKitXTB:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RDKitXTB)
include(CMakePackageConfigHelpers)
configure_package_config_file(cpp/rdXTBTools/RDKitXTBConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/RDKitXTBConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RDKitXTB)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/RDKitXTBConfigVersion.cmake"
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/RDKitXTBConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/RDKitXTBConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RDKitXTB)

if(BUILD_TESTING)
  add_executable(rdxtb_cpp_test cpp/rdXTBTools/testXTB.cpp)
  target_link_libraries(rdxtb_cpp_test PRIVATE RDKitXTB)
  add_test(NAME rdXTBTools/core COMMAND rdxtb_cpp_test)
  set_tests_properties(rdXTBTools/core PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
