# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "C:/Program Files/CMake/bin/cmake.exe"

# The command to remove a file.
RM = "C:/Program Files/CMake/bin/cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:/Users/Lenovo/Desktop/thermal/workstation/test0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:/Users/Lenovo/Desktop/thermal/workstation/test0/build

# Include any dependencies generated for this target.
include src/CMakeFiles/solid.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/solid.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/solid.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/solid.dir/flags.make

src/CMakeFiles/solid.dir/program_main.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/program_main.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/program_main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/CMakeFiles/solid.dir/program_main.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/program_main.f90 -o CMakeFiles/solid.dir/program_main.f90.obj

src/CMakeFiles/solid.dir/program_main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/program_main.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/program_main.f90 > CMakeFiles/solid.dir/program_main.f90.i

src/CMakeFiles/solid.dir/program_main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/program_main.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/program_main.f90 -o CMakeFiles/solid.dir/program_main.f90.s

src/CMakeFiles/solid.dir/dataprocess.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/dataprocess.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/dataprocess.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object src/CMakeFiles/solid.dir/dataprocess.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/dataprocess.f90 -o CMakeFiles/solid.dir/dataprocess.f90.obj

src/CMakeFiles/solid.dir/dataprocess.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/dataprocess.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/dataprocess.f90 > CMakeFiles/solid.dir/dataprocess.f90.i

src/CMakeFiles/solid.dir/dataprocess.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/dataprocess.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/dataprocess.f90 -o CMakeFiles/solid.dir/dataprocess.f90.s

src/CMakeFiles/solid.dir/driverthinput_test1.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/driverthinput_test1.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverthinput_test1.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object src/CMakeFiles/solid.dir/driverthinput_test1.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverthinput_test1.f90 -o CMakeFiles/solid.dir/driverthinput_test1.f90.obj

src/CMakeFiles/solid.dir/driverthinput_test1.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/driverthinput_test1.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverthinput_test1.f90 > CMakeFiles/solid.dir/driverthinput_test1.f90.i

src/CMakeFiles/solid.dir/driverthinput_test1.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/driverthinput_test1.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverthinput_test1.f90 -o CMakeFiles/solid.dir/driverthinput_test1.f90.s

src/CMakeFiles/solid.dir/driverTSoutput.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/driverTSoutput.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverTSoutput.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object src/CMakeFiles/solid.dir/driverTSoutput.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverTSoutput.f90 -o CMakeFiles/solid.dir/driverTSoutput.f90.obj

src/CMakeFiles/solid.dir/driverTSoutput.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/driverTSoutput.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverTSoutput.f90 > CMakeFiles/solid.dir/driverTSoutput.f90.i

src/CMakeFiles/solid.dir/driverTSoutput.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/driverTSoutput.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverTSoutput.f90 -o CMakeFiles/solid.dir/driverTSoutput.f90.s

src/CMakeFiles/solid.dir/globalTSconstants.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/globalTSconstants.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSconstants.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object src/CMakeFiles/solid.dir/globalTSconstants.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSconstants.f90 -o CMakeFiles/solid.dir/globalTSconstants.f90.obj

src/CMakeFiles/solid.dir/globalTSconstants.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/globalTSconstants.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSconstants.f90 > CMakeFiles/solid.dir/globalTSconstants.f90.i

src/CMakeFiles/solid.dir/globalTSconstants.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/globalTSconstants.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSconstants.f90 -o CMakeFiles/solid.dir/globalTSconstants.f90.s

src/CMakeFiles/solid.dir/globalTSvariables.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/globalTSvariables.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSvariables.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object src/CMakeFiles/solid.dir/globalTSvariables.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSvariables.f90 -o CMakeFiles/solid.dir/globalTSvariables.f90.obj

src/CMakeFiles/solid.dir/globalTSvariables.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/globalTSvariables.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSvariables.f90 > CMakeFiles/solid.dir/globalTSvariables.f90.i

src/CMakeFiles/solid.dir/globalTSvariables.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/globalTSvariables.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/globalTSvariables.f90 -o CMakeFiles/solid.dir/globalTSvariables.f90.s

src/CMakeFiles/solid.dir/driverSteadyTSsolver.f90.obj: src/CMakeFiles/solid.dir/flags.make
src/CMakeFiles/solid.dir/driverSteadyTSsolver.f90.obj: C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverSteadyTSsolver.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building Fortran object src/CMakeFiles/solid.dir/driverSteadyTSsolver.f90.obj"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverSteadyTSsolver.f90 -o CMakeFiles/solid.dir/driverSteadyTSsolver.f90.obj

src/CMakeFiles/solid.dir/driverSteadyTSsolver.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing Fortran source to CMakeFiles/solid.dir/driverSteadyTSsolver.f90.i"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverSteadyTSsolver.f90 > CMakeFiles/solid.dir/driverSteadyTSsolver.f90.i

src/CMakeFiles/solid.dir/driverSteadyTSsolver.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling Fortran source to assembly CMakeFiles/solid.dir/driverSteadyTSsolver.f90.s"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:/Users/Lenovo/Desktop/thermal/workstation/test0/src/driverSteadyTSsolver.f90 -o CMakeFiles/solid.dir/driverSteadyTSsolver.f90.s

# Object files for target solid
solid_OBJECTS = \
"CMakeFiles/solid.dir/program_main.f90.obj" \
"CMakeFiles/solid.dir/dataprocess.f90.obj" \
"CMakeFiles/solid.dir/driverthinput_test1.f90.obj" \
"CMakeFiles/solid.dir/driverTSoutput.f90.obj" \
"CMakeFiles/solid.dir/globalTSconstants.f90.obj" \
"CMakeFiles/solid.dir/globalTSvariables.f90.obj" \
"CMakeFiles/solid.dir/driverSteadyTSsolver.f90.obj"

# External object files for target solid
solid_EXTERNAL_OBJECTS =

src/solid.exe: src/CMakeFiles/solid.dir/program_main.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/dataprocess.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/driverthinput_test1.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/driverTSoutput.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/globalTSconstants.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/globalTSvariables.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/driverSteadyTSsolver.f90.obj
src/solid.exe: src/CMakeFiles/solid.dir/build.make
src/solid.exe: src/CMakeFiles/solid.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking Fortran executable solid.exe"
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && "C:/Program Files/CMake/bin/cmake.exe" -E rm -f CMakeFiles/solid.dir/objects.a
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/ar.exe qc CMakeFiles/solid.dir/objects.a @CMakeFiles/solid.dir/objects1.rsp
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && C:/MinGW/bin/gfortran.exe -g -Wl,--whole-archive CMakeFiles/solid.dir/objects.a -Wl,--no-whole-archive -o solid.exe -Wl,--out-implib,libsolid.dll.a -Wl,--major-image-version,0,--minor-image-version,0 

# Rule to build all files generated by this target.
src/CMakeFiles/solid.dir/build: src/solid.exe
.PHONY : src/CMakeFiles/solid.dir/build

src/CMakeFiles/solid.dir/clean:
	cd C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src && $(CMAKE_COMMAND) -P CMakeFiles/solid.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/solid.dir/clean

src/CMakeFiles/solid.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" C:/Users/Lenovo/Desktop/thermal/workstation/test0 C:/Users/Lenovo/Desktop/thermal/workstation/test0/src C:/Users/Lenovo/Desktop/thermal/workstation/test0/build C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src C:/Users/Lenovo/Desktop/thermal/workstation/test0/build/src/CMakeFiles/solid.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : src/CMakeFiles/solid.dir/depend

