# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/daq/sw_daq_tofhir_v1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/daq/sw_daq_tofhir_v1/build

# Include any dependencies generated for this target.
include CMakeFiles/daqd.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/daqd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/daqd.dir/flags.make

CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o: CMakeFiles/daqd.dir/flags.make
CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o: ../src/daqd/daqd.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/daqd/daqd.cpp

CMakeFiles/daqd.dir/src/daqd/daqd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqd.dir/src/daqd/daqd.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/daqd/daqd.cpp > CMakeFiles/daqd.dir/src/daqd/daqd.cpp.i

CMakeFiles/daqd.dir/src/daqd/daqd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqd.dir/src/daqd/daqd.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/daqd/daqd.cpp -o CMakeFiles/daqd.dir/src/daqd/daqd.cpp.s

CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.requires:
.PHONY : CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.requires

CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.provides: CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.requires
	$(MAKE) -f CMakeFiles/daqd.dir/build.make CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.provides.build
.PHONY : CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.provides

CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.provides.build: CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o

CMakeFiles/daqd.dir/src/daqd/Client.cpp.o: CMakeFiles/daqd.dir/flags.make
CMakeFiles/daqd.dir/src/daqd/Client.cpp.o: ../src/daqd/Client.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/daqd.dir/src/daqd/Client.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/daqd.dir/src/daqd/Client.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/daqd/Client.cpp

CMakeFiles/daqd.dir/src/daqd/Client.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqd.dir/src/daqd/Client.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/daqd/Client.cpp > CMakeFiles/daqd.dir/src/daqd/Client.cpp.i

CMakeFiles/daqd.dir/src/daqd/Client.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqd.dir/src/daqd/Client.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/daqd/Client.cpp -o CMakeFiles/daqd.dir/src/daqd/Client.cpp.s

CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.requires:
.PHONY : CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.requires

CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.provides: CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.requires
	$(MAKE) -f CMakeFiles/daqd.dir/build.make CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.provides.build
.PHONY : CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.provides

CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.provides.build: CMakeFiles/daqd.dir/src/daqd/Client.cpp.o

CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o: CMakeFiles/daqd.dir/flags.make
CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o: ../src/daqd/FrameServer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/daqd/FrameServer.cpp

CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/daqd/FrameServer.cpp > CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.i

CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/daqd/FrameServer.cpp -o CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.s

CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.requires:
.PHONY : CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.requires

CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.provides: CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.requires
	$(MAKE) -f CMakeFiles/daqd.dir/build.make CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.provides.build
.PHONY : CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.provides

CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.provides.build: CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o

CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o: CMakeFiles/daqd.dir/flags.make
CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o: ../src/daqd/UDPFrameServer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/daqd/UDPFrameServer.cpp

CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/daqd/UDPFrameServer.cpp > CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.i

CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/daqd/UDPFrameServer.cpp -o CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.s

CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.requires:
.PHONY : CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.requires

CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.provides: CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.requires
	$(MAKE) -f CMakeFiles/daqd.dir/build.make CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.provides.build
.PHONY : CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.provides

CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.provides.build: CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o

CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o: CMakeFiles/daqd.dir/flags.make
CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o: ../src/raw_data/shm_raw.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/raw_data/shm_raw.cpp

CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/raw_data/shm_raw.cpp > CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.i

CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/raw_data/shm_raw.cpp -o CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.s

CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.requires:
.PHONY : CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.requires

CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.provides: CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.requires
	$(MAKE) -f CMakeFiles/daqd.dir/build.make CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.provides.build
.PHONY : CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.provides

CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.provides.build: CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o

# Object files for target daqd
daqd_OBJECTS = \
"CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o" \
"CMakeFiles/daqd.dir/src/daqd/Client.cpp.o" \
"CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o" \
"CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o" \
"CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o"

# External object files for target daqd
daqd_EXTERNAL_OBJECTS =

daqd: CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o
daqd: CMakeFiles/daqd.dir/src/daqd/Client.cpp.o
daqd: CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o
daqd: CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o
daqd: CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o
daqd: CMakeFiles/daqd.dir/build.make
daqd: CMakeFiles/daqd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable daqd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/daqd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/daqd.dir/build: daqd
.PHONY : CMakeFiles/daqd.dir/build

CMakeFiles/daqd.dir/requires: CMakeFiles/daqd.dir/src/daqd/daqd.cpp.o.requires
CMakeFiles/daqd.dir/requires: CMakeFiles/daqd.dir/src/daqd/Client.cpp.o.requires
CMakeFiles/daqd.dir/requires: CMakeFiles/daqd.dir/src/daqd/FrameServer.cpp.o.requires
CMakeFiles/daqd.dir/requires: CMakeFiles/daqd.dir/src/daqd/UDPFrameServer.cpp.o.requires
CMakeFiles/daqd.dir/requires: CMakeFiles/daqd.dir/src/raw_data/shm_raw.cpp.o.requires
.PHONY : CMakeFiles/daqd.dir/requires

CMakeFiles/daqd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/daqd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/daqd.dir/clean

CMakeFiles/daqd.dir/depend:
	cd /home/daq/sw_daq_tofhir_v1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/daq/sw_daq_tofhir_v1 /home/daq/sw_daq_tofhir_v1 /home/daq/sw_daq_tofhir_v1/build /home/daq/sw_daq_tofhir_v1/build /home/daq/sw_daq_tofhir_v1/build/CMakeFiles/daqd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/daqd.dir/depend

