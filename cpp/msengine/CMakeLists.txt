cmake_minimum_required(VERSION 3.16)
project(msengine CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Test builds only (the release library is built with the defaults):
#   -DMSENGINE_SANITIZE=address,undefined   AddressSanitizer + UBSan (+ float-cast-overflow), -O1 -g
#   -DMSENGINE_SANITIZE=thread              ThreadSanitizer (use -DMSENGINE_OPENMP=OFF: libgomp is not instrumented)
#   -DMSENGINE_OPENMP=OFF                   no OpenMP (the pragmas are ignored, one thread)
option(MSENGINE_OPENMP "build with OpenMP" ON)
set(MSENGINE_SANITIZE "" CACHE STRING "sanitizers of a test build (address,undefined or thread)")

if(MSENGINE_OPENMP)
  find_package(OpenMP)
endif()
set(SOURCES
  src/common.cpp
  src/spectrum.cpp
  src/mzml.cpp
  src/bruker.cpp
  src/file.cpp
  src/peaks.cpp
  src/maxent.cpp
  src/shimadzu.cpp
  src/arrays.cpp
  src/lctools.cpp
  src/calib.cpp
  src/formula.cpp
  src/unidec.cpp
  src/isodec.cpp
  src/udengine/udengine.cpp
  src/udengine/udscore.cpp
  src/pda.cpp
  src/compare.cpp
  src/adducts.cpp
  src/shifts.cpp
  src/vendor.cpp
  third_party/miniz.c
  third_party/sqlite3.c)
set(MS_DEFS MSENGINE_BUILD SQLITE_OMIT_LOAD_EXTENSION SQLITE_THREADSAFE=1 SQLITE_DEFAULT_MEMSTATUS=0)

# the sources are compiled once and linked into the library and into the reader test
add_library(msengine_obj OBJECT ${SOURCES})
set_target_properties(msengine_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(msengine_obj PUBLIC include PRIVATE src third_party)
target_compile_definitions(msengine_obj PRIVATE ${MS_DEFS})
set_source_files_properties(third_party/sqlite3.c PROPERTIES COMPILE_OPTIONS "-w")
# the UniDec engine: single precision arithmetic as written (no fused multiply add on any target), so that
# every build gives the same results (tests/check_udengine.py)
set_source_files_properties(src/udengine/udengine.cpp src/udengine/udscore.cpp PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
# the mass shift finder: the arithmetic of its Python reference (ms_shifts.py), no fused multiply add
set_source_files_properties(src/shifts.cpp PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
target_compile_options(msengine_obj PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O3 -Wall -Wextra>)
if(NOT WIN32)
  # only the ms_* functions are exported (MS_API): otherwise the vendor libraries loaded later
  # (timsdata, baf2sql) bind their own inline C++ functions to this library's copies
  target_compile_options(msengine_obj PRIVATE -fvisibility=hidden $<$<COMPILE_LANGUAGE:CXX>:-fvisibility-inlines-hidden>)
endif()

add_library(msengine SHARED $<TARGET_OBJECTS:msengine_obj>)
target_include_directories(msengine PUBLIC include)
set_target_properties(msengine PROPERTIES PREFIX "" OUTPUT_NAME "msengine")
if(NOT WIN32 AND NOT APPLE)
  target_link_options(msengine PRIVATE -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/cmake/msengine.map)
  set_property(TARGET msengine APPEND PROPERTY LINK_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/msengine.map)
endif()

# test programs: only when the tests folder is present (a copy without tests builds the library alone)
set(MS_TESTS OFF)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/ms_test.cpp)
  set(MS_TESTS ON)
endif()
set(MS_TARGETS msengine)
set(MS_OBJ_TESTS)
if(MS_TESTS)
# command line tests (compare with the Python reference)
add_executable(ms_test tests/ms_test.cpp)
target_link_libraries(ms_test msengine)

# reader level test binary (uses ms::Reader directly, not exported from the DLL)
add_executable(ms_test_e3 tests/ms_test.cpp $<TARGET_OBJECTS:msengine_obj>)
target_include_directories(ms_test_e3 PRIVATE include src third_party)
target_compile_definitions(ms_test_e3 PRIVATE ${MS_DEFS} MS_TEST_E3)

# robustness tests of the C API (tests/README.md): bad arguments, damaged files, threads, cancellation
add_executable(robust_test tests/robust_test.cpp)
target_link_libraries(robust_test msengine)

# the UniDec engine of the library on a conf.dat, for the comparison with unidec.exe (tests/check_udengine.py)
add_executable(ud_engine_test tests/ud_engine_test.cpp $<TARGET_OBJECTS:msengine_obj>)
target_include_directories(ud_engine_test PRIVATE include src third_party)
target_compile_definitions(ud_engine_test PRIVATE ${MS_DEFS})

# threads, cancellation and the stops of the UniDec engine in the library
add_executable(udengine_test tests/udengine_test.cpp $<TARGET_OBJECTS:msengine_obj>)
target_include_directories(udengine_test PRIVATE include src third_party)
target_compile_definitions(udengine_test PRIVATE ${MS_DEFS})

list(APPEND MS_TARGETS ms_test ms_test_e3 robust_test ud_engine_test udengine_test)
set(MS_OBJ_TESTS ms_test_e3 ud_engine_test udengine_test)
endif()

if(OpenMP_CXX_FOUND)
  if(MINGW)
    # FindOpenMP picks the import library of libgomp-1.dll; -fopenmp with the
    # static link flags of the toolchain pulls libgomp.a into the DLL instead
    target_compile_options(msengine_obj PRIVATE -fopenmp)
    foreach(t msengine ${MS_OBJ_TESTS})
      target_link_options(${t} PRIVATE -fopenmp)
    endforeach()
  else()
    target_link_libraries(msengine_obj PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(msengine PUBLIC OpenMP::OpenMP_CXX)
    foreach(t ${MS_OBJ_TESTS})
      target_link_libraries(${t} OpenMP::OpenMP_CXX)
    endforeach()
  endif()
endif()
if(NOT WIN32)
  target_link_libraries(msengine PUBLIC dl)
  if(MS_TESTS)
    find_package(Threads REQUIRED)
    foreach(t ${MS_OBJ_TESTS})
      target_link_libraries(${t} dl Threads::Threads)
    endforeach()
    target_link_libraries(robust_test Threads::Threads)
  endif()
endif()

if(MSENGINE_SANITIZE)
  set(SAN_FLAGS -fsanitize=${MSENGINE_SANITIZE} -fno-omit-frame-pointer -g -O1)
  if(MSENGINE_SANITIZE MATCHES "undefined")
    list(APPEND SAN_FLAGS -fsanitize=float-cast-overflow)
  endif()
  target_compile_options(msengine_obj PRIVATE ${SAN_FLAGS})
  foreach(t ${MS_TARGETS})
    target_compile_options(${t} PRIVATE ${SAN_FLAGS})
    target_link_options(${t} PRIVATE -fsanitize=${MSENGINE_SANITIZE})
  endforeach()
  # SQLite (third party, 9 MB of C) is not instrumented: it builds slowly and is fuzzed upstream
  set_source_files_properties(third_party/sqlite3.c PROPERTIES COMPILE_OPTIONS "-w;-fno-sanitize=all;-O2")
endif()
