add_library(
    psi
    OBJECT
    psi.cpp
)

add_library(
  psi_initializer
  OBJECT
  psi_initializer.cpp
  psi_initializer_atomic.cpp
  psi_initializer_random.cpp
  psi_initializer_nao.cpp
)

if(ENABLE_COVERAGE)
  add_coverage(psi)
endif()

if (BUILD_TESTING)
  add_subdirectory(kernels/test)
  add_subdirectory(test)
endif()