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
include CMakeFiles/common.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/common.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/common.dir/flags.make

CMakeFiles/common.dir/src/base/ThreadPool.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/ThreadPool.cpp.o: ../src/base/ThreadPool.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/ThreadPool.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/ThreadPool.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/ThreadPool.cpp

CMakeFiles/common.dir/src/base/ThreadPool.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/ThreadPool.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/ThreadPool.cpp > CMakeFiles/common.dir/src/base/ThreadPool.cpp.i

CMakeFiles/common.dir/src/base/ThreadPool.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/ThreadPool.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/ThreadPool.cpp -o CMakeFiles/common.dir/src/base/ThreadPool.cpp.s

CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.requires

CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.provides: CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.provides

CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.provides.build: CMakeFiles/common.dir/src/base/ThreadPool.cpp.o

CMakeFiles/common.dir/src/base/SystemConfig.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/SystemConfig.cpp.o: ../src/base/SystemConfig.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/SystemConfig.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/SystemConfig.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/SystemConfig.cpp

CMakeFiles/common.dir/src/base/SystemConfig.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/SystemConfig.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/SystemConfig.cpp > CMakeFiles/common.dir/src/base/SystemConfig.cpp.i

CMakeFiles/common.dir/src/base/SystemConfig.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/SystemConfig.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/SystemConfig.cpp -o CMakeFiles/common.dir/src/base/SystemConfig.cpp.s

CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.requires

CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.provides: CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.provides

CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.provides.build: CMakeFiles/common.dir/src/base/SystemConfig.cpp.o

CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o: ../src/raw_data/RawReader.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/raw_data/RawReader.cpp

CMakeFiles/common.dir/src/raw_data/RawReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/raw_data/RawReader.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/raw_data/RawReader.cpp > CMakeFiles/common.dir/src/raw_data/RawReader.cpp.i

CMakeFiles/common.dir/src/raw_data/RawReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/raw_data/RawReader.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/raw_data/RawReader.cpp -o CMakeFiles/common.dir/src/raw_data/RawReader.cpp.s

CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.requires

CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.provides: CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.provides

CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.provides.build: CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o

CMakeFiles/common.dir/src/base/Instrumentation.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/Instrumentation.cpp.o: ../src/base/Instrumentation.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/Instrumentation.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/Instrumentation.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/Instrumentation.cpp

CMakeFiles/common.dir/src/base/Instrumentation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/Instrumentation.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/Instrumentation.cpp > CMakeFiles/common.dir/src/base/Instrumentation.cpp.i

CMakeFiles/common.dir/src/base/Instrumentation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/Instrumentation.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/Instrumentation.cpp -o CMakeFiles/common.dir/src/base/Instrumentation.cpp.s

CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.requires

CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.provides: CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.provides

CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.provides.build: CMakeFiles/common.dir/src/base/Instrumentation.cpp.o

CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o: ../src/base/CoarseSorter.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/CoarseSorter.cpp

CMakeFiles/common.dir/src/base/CoarseSorter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/CoarseSorter.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/CoarseSorter.cpp > CMakeFiles/common.dir/src/base/CoarseSorter.cpp.i

CMakeFiles/common.dir/src/base/CoarseSorter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/CoarseSorter.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/CoarseSorter.cpp -o CMakeFiles/common.dir/src/base/CoarseSorter.cpp.s

CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.requires

CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.provides: CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.provides

CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.provides.build: CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o

CMakeFiles/common.dir/src/base/ProcessHit.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/ProcessHit.cpp.o: ../src/base/ProcessHit.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/ProcessHit.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/ProcessHit.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/ProcessHit.cpp

CMakeFiles/common.dir/src/base/ProcessHit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/ProcessHit.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/ProcessHit.cpp > CMakeFiles/common.dir/src/base/ProcessHit.cpp.i

CMakeFiles/common.dir/src/base/ProcessHit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/ProcessHit.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/ProcessHit.cpp -o CMakeFiles/common.dir/src/base/ProcessHit.cpp.s

CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.requires

CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.provides: CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.provides

CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.provides.build: CMakeFiles/common.dir/src/base/ProcessHit.cpp.o

CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o: ../src/base/SimpleGrouper.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/SimpleGrouper.cpp

CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/SimpleGrouper.cpp > CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.i

CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/SimpleGrouper.cpp -o CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.s

CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.requires

CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.provides: CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.provides

CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.provides.build: CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o

CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o: CMakeFiles/common.dir/flags.make
CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o: ../src/base/CoincidenceGrouper.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/daq/sw_daq_tofhir_v1/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o -c /home/daq/sw_daq_tofhir_v1/src/base/CoincidenceGrouper.cpp

CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/daq/sw_daq_tofhir_v1/src/base/CoincidenceGrouper.cpp > CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.i

CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/daq/sw_daq_tofhir_v1/src/base/CoincidenceGrouper.cpp -o CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.s

CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.requires:
.PHONY : CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.requires

CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.provides: CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.requires
	$(MAKE) -f CMakeFiles/common.dir/build.make CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.provides.build
.PHONY : CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.provides

CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.provides.build: CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o

# Object files for target common
common_OBJECTS = \
"CMakeFiles/common.dir/src/base/ThreadPool.cpp.o" \
"CMakeFiles/common.dir/src/base/SystemConfig.cpp.o" \
"CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o" \
"CMakeFiles/common.dir/src/base/Instrumentation.cpp.o" \
"CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o" \
"CMakeFiles/common.dir/src/base/ProcessHit.cpp.o" \
"CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o" \
"CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o"

# External object files for target common
common_EXTERNAL_OBJECTS =

libcommon.a: CMakeFiles/common.dir/src/base/ThreadPool.cpp.o
libcommon.a: CMakeFiles/common.dir/src/base/SystemConfig.cpp.o
libcommon.a: CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o
libcommon.a: CMakeFiles/common.dir/src/base/Instrumentation.cpp.o
libcommon.a: CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o
libcommon.a: CMakeFiles/common.dir/src/base/ProcessHit.cpp.o
libcommon.a: CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o
libcommon.a: CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o
libcommon.a: CMakeFiles/common.dir/build.make
libcommon.a: CMakeFiles/common.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libcommon.a"
	$(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/common.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/common.dir/build: libcommon.a
.PHONY : CMakeFiles/common.dir/build

CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/ThreadPool.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/SystemConfig.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/raw_data/RawReader.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/Instrumentation.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/CoarseSorter.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/ProcessHit.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/SimpleGrouper.cpp.o.requires
CMakeFiles/common.dir/requires: CMakeFiles/common.dir/src/base/CoincidenceGrouper.cpp.o.requires
.PHONY : CMakeFiles/common.dir/requires

CMakeFiles/common.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/common.dir/cmake_clean.cmake
.PHONY : CMakeFiles/common.dir/clean

CMakeFiles/common.dir/depend:
	cd /home/daq/sw_daq_tofhir_v1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/daq/sw_daq_tofhir_v1 /home/daq/sw_daq_tofhir_v1 /home/daq/sw_daq_tofhir_v1/build /home/daq/sw_daq_tofhir_v1/build /home/daq/sw_daq_tofhir_v1/build/CMakeFiles/common.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/common.dir/depend

