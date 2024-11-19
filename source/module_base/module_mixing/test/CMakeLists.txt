# the following test should be run with only one omp thread
# use googletest

add_executable(test_mixing mixing_test.cpp)
target_link_libraries(test_mixing parameter base device ${math_libs}  Threads::Threads
GTest::gtest_main GTest::gmock_main)


add_test(
  NAME test_mixing
  COMMAND OMP_NUM_THREADS=1 mpirun -np 1 test_mixing
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# AddTest(
#   TARGET test_mixing
#   LIBS parameter base device ${math_libs} 
#   SOURCES mixing_test.cpp
# )