include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# target 
file( GLOB MPC_SRCS *.cpp *.hpp)
add_executable(mpc ${MPC_SRCS})           
target_link_libraries(mpc ${ARMADILLO_LIBRARIES}) 