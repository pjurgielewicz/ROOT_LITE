############################################################################
# CMakeLists.txt file for building ROOT core/cont package
############################################################################

if(testing)
  add_subdirectory(test)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/res ${CMAKE_CURRENT_SOURCE_DIR}/../foundation/res)

ROOT_GLOB_HEADERS(headers ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)
set(Cont_dict_headers ${headers} PARENT_SCOPE)

ROOT_OBJECT_LIBRARY(Cont *.cxx)

ROOT_INSTALL_HEADERS()
