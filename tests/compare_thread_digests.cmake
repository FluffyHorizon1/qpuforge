# Runs test_determinism --digest under OMP_NUM_THREADS = 1, 2, 4, 8 (separate processes, so the environment variable
# itself is what changes) and fails unless all four outputs are identical.
#   cmake -DEXE=<path to test_determinism> -P compare_thread_digests.cmake
set(_ref "")
foreach(n 1 2 4 8)
  execute_process(COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${n} ${EXE} --digest
                  OUTPUT_VARIABLE _out RESULT_VARIABLE _rc OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "test_determinism --digest failed at OMP_NUM_THREADS=${n} (exit ${_rc})")
  endif()
  if(n EQUAL 1)
    set(_ref "${_out}")
    message(STATUS "OMP_NUM_THREADS=1: ${_out}")
  elseif(NOT _out STREQUAL _ref)
    message(FATAL_ERROR "OMP_NUM_THREADS=${n} differs from OMP_NUM_THREADS=1:\n  1: ${_ref}\n  ${n}: ${_out}")
  else()
    message(STATUS "OMP_NUM_THREADS=${n}: identical")
  endif()
endforeach()
