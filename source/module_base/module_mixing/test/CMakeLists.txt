# it is nonsence to disable MPI for so-called simplicitly of unittest
AddTest(
  TARGET test_mixing
  LIBS parameter base device ${math_libs} 
  SOURCES mixing_test.cpp
)