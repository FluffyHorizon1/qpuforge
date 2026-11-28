# The grid solve of flux_coupling against tools/flux_check.py, a separate numpy implementation of the same
# discretisation (Jacobi preconditioner; the two share no code): AF-17 from the development cells, the data qubit Q8 and
# its line, on 1 um cells. layout_gates dumps the window and its result; flux_check.py solves the dump again and fails
# when M or M_cm differs by more than 0.001 pH.
# Usage: cmake -DEXE=<layout_gates> -DPYTHON=<python3 with numpy> -DTOOL=<tools/flux_check.py> -DCELLS=<cells.json>
#              -DOUT=<window.json> -P check_flux_against_numpy.cmake
foreach(v EXE PYTHON TOOL CELLS OUT)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "check_flux_against_numpy.cmake: -D${v}=... is required")
  endif()
endforeach()
execute_process(COMMAND "${EXE}" flux alpha1-v11 "${CELLS}" Q8 FQ8 --cell 1 --dump "${OUT}" RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "layout_gates flux: exit status ${rc}\n${out}\n${err}")
endif()
message(STATUS "${out}")
execute_process(COMMAND "${PYTHON}" -I "${TOOL}" "${OUT}" --h 1 --max-diff 0.001 RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}")
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "tools/flux_check.py: exit status ${rc}\n${err}")
endif()
