include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# target 
file( GLOB CLQR_SRCS *.cpp *.hpp)
add_executable(clqr ${CLQR_SRCS})           
target_link_libraries(clqr ${ARMADILLO_LIBRARIES}) 