# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lur/Desktop/inviwo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lur/Desktop/inviwo

# Include any dependencies generated for this target.
include modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/depend.make

# Include the progress variables for this target.
include modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/progress.make

# Include the compile flags for this target's objects.
include modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/flags.make

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.o: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/flags.make
modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.o: modules/labmarchingsquares/marchingsquares.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lur/Desktop/inviwo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.o"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.o -c /home/lur/Desktop/inviwo/modules/labmarchingsquares/marchingsquares.cpp

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.i"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lur/Desktop/inviwo/modules/labmarchingsquares/marchingsquares.cpp > CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.i

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.s"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lur/Desktop/inviwo/modules/labmarchingsquares/marchingsquares.cpp -o CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.s

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.o: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/flags.make
modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.o: modules/labmarchingsquares/utils/amirameshvolumereader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lur/Desktop/inviwo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.o"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.o -c /home/lur/Desktop/inviwo/modules/labmarchingsquares/utils/amirameshvolumereader.cpp

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.i"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lur/Desktop/inviwo/modules/labmarchingsquares/utils/amirameshvolumereader.cpp > CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.i

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.s"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lur/Desktop/inviwo/modules/labmarchingsquares/utils/amirameshvolumereader.cpp -o CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.s

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.o: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/flags.make
modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.o: modules/labmarchingsquares/utils/setminmaxdatamap.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lur/Desktop/inviwo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.o"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.o -c /home/lur/Desktop/inviwo/modules/labmarchingsquares/utils/setminmaxdatamap.cpp

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.i"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lur/Desktop/inviwo/modules/labmarchingsquares/utils/setminmaxdatamap.cpp > CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.i

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.s"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lur/Desktop/inviwo/modules/labmarchingsquares/utils/setminmaxdatamap.cpp -o CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.s

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.o: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/flags.make
modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.o: modules/labmarchingsquares/labmarchingsquaresmodule.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lur/Desktop/inviwo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.o"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.o -c /home/lur/Desktop/inviwo/modules/labmarchingsquares/labmarchingsquaresmodule.cpp

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.i"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lur/Desktop/inviwo/modules/labmarchingsquares/labmarchingsquaresmodule.cpp > CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.i

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.s"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lur/Desktop/inviwo/modules/labmarchingsquares/labmarchingsquaresmodule.cpp -o CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.s

# Object files for target inviwo-module-labmarchingsquares
inviwo__module__labmarchingsquares_OBJECTS = \
"CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.o" \
"CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.o" \
"CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.o" \
"CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.o"

# External object files for target inviwo-module-labmarchingsquares
inviwo__module__labmarchingsquares_EXTERNAL_OBJECTS =

lib/libinviwo-module-labmarchingsquares.so: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/marchingsquares.cpp.o
lib/libinviwo-module-labmarchingsquares.so: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/amirameshvolumereader.cpp.o
lib/libinviwo-module-labmarchingsquares.so: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/utils/setminmaxdatamap.cpp.o
lib/libinviwo-module-labmarchingsquares.so: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/labmarchingsquaresmodule.cpp.o
lib/libinviwo-module-labmarchingsquares.so: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/build.make
lib/libinviwo-module-labmarchingsquares.so: lib/libinviwo-core.so
lib/libinviwo-module-labmarchingsquares.so: lib/libticpp.so
lib/libinviwo-module-labmarchingsquares.so: lib/libsigar.so
lib/libinviwo-module-labmarchingsquares.so: modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lur/Desktop/inviwo/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX shared library ../../lib/libinviwo-module-labmarchingsquares.so"
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/inviwo-module-labmarchingsquares.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/build: lib/libinviwo-module-labmarchingsquares.so

.PHONY : modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/build

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/clean:
	cd /home/lur/Desktop/inviwo/modules/labmarchingsquares && $(CMAKE_COMMAND) -P CMakeFiles/inviwo-module-labmarchingsquares.dir/cmake_clean.cmake
.PHONY : modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/clean

modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/depend:
	cd /home/lur/Desktop/inviwo && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lur/Desktop/inviwo /home/lur/Desktop/inviwo/modules/labmarchingsquares /home/lur/Desktop/inviwo /home/lur/Desktop/inviwo/modules/labmarchingsquares /home/lur/Desktop/inviwo/modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : modules/labmarchingsquares/CMakeFiles/inviwo-module-labmarchingsquares.dir/depend

