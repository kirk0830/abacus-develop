AddTest(
    TARGET psi_UT
    LIBS ${math_libs} base device
    SOURCES 
        psi_test.cpp
        ../psi.cpp 
)

if(ENABLE_LCAO)
AddTest(
    TARGET psi_initializer_integrate_test
    SOURCES
        integrated_initializer_test.cpp
)
endif()
