SET(SOURCE_FILES
    test_solver.f90
	../src/dataprocess.f90
	../src/driverthinput_test1.f90
	../src/driverTSoutput.f90
	../src/globalTSconstants.f90
	../src/globalTSvariables.f90
	../src/driverSteadyTSsolver.f90
    )

ADD_EXECUTABLE(test_solver ${SOURCE_FILES})

# SET(SOURCE_FILES

#     )

# ADD_EXECUTABLE(test_cuboid ${SOURCE_FILES})

# SET(SOURCE_FILES

#     )

# ADD_EXECUTABLE(test_triangular_prism ${SOURCE_FILES})

ADD_TEST(NAME unitTest_solid COMMAND test_solid)