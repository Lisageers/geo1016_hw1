# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug

# Include any dependencies generated for this target.
include easy3d/core/CMakeFiles/easy3d_core.dir/depend.make

# Include the progress variables for this target.
include easy3d/core/CMakeFiles/easy3d_core.dir/progress.make

# Include the compile flags for this target's objects.
include easy3d/core/CMakeFiles/easy3d_core.dir/flags.make

easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.o: easy3d/core/CMakeFiles/easy3d_core.dir/flags.make
easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.o: ../easy3d/core/graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.o"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/easy3d_core.dir/graph.cpp.o -c /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/graph.cpp

easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/graph.cpp.i"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/graph.cpp > CMakeFiles/easy3d_core.dir/graph.cpp.i

easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/graph.cpp.s"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/graph.cpp -o CMakeFiles/easy3d_core.dir/graph.cpp.s

easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.o: easy3d/core/CMakeFiles/easy3d_core.dir/flags.make
easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.o: ../easy3d/core/kdtree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.o"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/easy3d_core.dir/kdtree.cpp.o -c /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/kdtree.cpp

easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/kdtree.cpp.i"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/kdtree.cpp > CMakeFiles/easy3d_core.dir/kdtree.cpp.i

easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/kdtree.cpp.s"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/kdtree.cpp -o CMakeFiles/easy3d_core.dir/kdtree.cpp.s

easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.o: easy3d/core/CMakeFiles/easy3d_core.dir/flags.make
easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.o: ../easy3d/core/point_cloud.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.o"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/easy3d_core.dir/point_cloud.cpp.o -c /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/point_cloud.cpp

easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/point_cloud.cpp.i"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/point_cloud.cpp > CMakeFiles/easy3d_core.dir/point_cloud.cpp.i

easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/point_cloud.cpp.s"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/point_cloud.cpp -o CMakeFiles/easy3d_core.dir/point_cloud.cpp.s

easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.o: easy3d/core/CMakeFiles/easy3d_core.dir/flags.make
easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.o: ../easy3d/core/surface_mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.o"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/easy3d_core.dir/surface_mesh.cpp.o -c /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/surface_mesh.cpp

easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/surface_mesh.cpp.i"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/surface_mesh.cpp > CMakeFiles/easy3d_core.dir/surface_mesh.cpp.i

easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/surface_mesh.cpp.s"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/surface_mesh.cpp -o CMakeFiles/easy3d_core.dir/surface_mesh.cpp.s

easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.o: easy3d/core/CMakeFiles/easy3d_core.dir/flags.make
easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.o: ../easy3d/core/manifold_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.o"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/easy3d_core.dir/manifold_builder.cpp.o -c /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/manifold_builder.cpp

easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/manifold_builder.cpp.i"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/manifold_builder.cpp > CMakeFiles/easy3d_core.dir/manifold_builder.cpp.i

easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/manifold_builder.cpp.s"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core/manifold_builder.cpp -o CMakeFiles/easy3d_core.dir/manifold_builder.cpp.s

# Object files for target easy3d_core
easy3d_core_OBJECTS = \
"CMakeFiles/easy3d_core.dir/graph.cpp.o" \
"CMakeFiles/easy3d_core.dir/kdtree.cpp.o" \
"CMakeFiles/easy3d_core.dir/point_cloud.cpp.o" \
"CMakeFiles/easy3d_core.dir/surface_mesh.cpp.o" \
"CMakeFiles/easy3d_core.dir/manifold_builder.cpp.o"

# External object files for target easy3d_core
easy3d_core_EXTERNAL_OBJECTS =

lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.o
lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.o
lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.o
lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.o
lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.o
lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/build.make
lib/libeasy3d_core.a: easy3d/core/CMakeFiles/easy3d_core.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library ../../lib/libeasy3d_core.a"
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && $(CMAKE_COMMAND) -P CMakeFiles/easy3d_core.dir/cmake_clean_target.cmake
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/easy3d_core.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
easy3d/core/CMakeFiles/easy3d_core.dir/build: lib/libeasy3d_core.a

.PHONY : easy3d/core/CMakeFiles/easy3d_core.dir/build

easy3d/core/CMakeFiles/easy3d_core.dir/clean:
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core && $(CMAKE_COMMAND) -P CMakeFiles/easy3d_core.dir/cmake_clean.cmake
.PHONY : easy3d/core/CMakeFiles/easy3d_core.dir/clean

easy3d/core/CMakeFiles/easy3d_core.dir/depend:
	cd /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/easy3d/core /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core /Users/Vince/Documents/Geo1016/A1_Calibration/A1_Calibration_Code/cmake-build-debug/easy3d/core/CMakeFiles/easy3d_core.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : easy3d/core/CMakeFiles/easy3d_core.dir/depend

