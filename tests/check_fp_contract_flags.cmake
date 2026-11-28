# q13: the FMA guard must share compile flags with the real targets.
# If test_fp_contract lists -ffp-contract=off privately, it cannot detect the flag
# being dropped from qpuforge / the review programs / test_qpuforge.
if(NOT SOURCE_DIR)
  message(FATAL_ERROR "SOURCE_DIR is not set")
endif()
file(READ "${SOURCE_DIR}/CMakeLists.txt" txt)
if(NOT txt MATCHES [[set\(QPUFORGE_COMPILE_OPTIONS]])
  message(FATAL_ERROR "CMakeLists.txt must define QPUFORGE_COMPILE_OPTIONS shared by every target")
endif()
if(NOT txt MATCHES [[QPUFORGE_COMPILE_OPTIONS[^\n]*-ffp-contract=off]])
  message(FATAL_ERROR "QPUFORGE_COMPILE_OPTIONS must include -ffp-contract=off")
endif()
if(txt MATCHES [[target_compile_options\(test_fp_contract[^\n]*-ffp-contract=]])
  message(FATAL_ERROR "test_fp_contract has its own -ffp-contract flag; use QPUFORGE_COMPILE_OPTIONS so the guard tracks the real targets")
endif()
if(NOT txt MATCHES [[target_compile_options\(\$\{name\} PRIVATE \$\{QPUFORGE_COMPILE_OPTIONS\}\)]]
    AND NOT txt MATCHES [[test_fp_contract:tests/test_fp_contract.cpp]])
  message(FATAL_ERROR "the target helper must use QPUFORGE_COMPILE_OPTIONS (or test_fp_contract must be built by that helper)")
endif()
if(NOT txt MATCHES [[test_fp_contract:tests/test_fp_contract.cpp]]
    AND NOT txt MATCHES [[target_compile_options\(test_fp_contract PRIVATE \$\{QPUFORGE_COMPILE_OPTIONS\}\)]])
  message(FATAL_ERROR "test_fp_contract is not built with QPUFORGE_COMPILE_OPTIONS")
endif()
# Every target that compiles something must be given the list by a line of its own. Each add_executable(<name> ...)
# and add_library(<name> ...) in CMakeLists.txt needs
#   target_compile_options(<name> PRIVATE ${QPUFORGE_COMPILE_OPTIONS})
# with the same <name>, whether it is written out (test_determinism, qpuforge_git_stamp) or a variable (${name} in the
# target helper). INTERFACE, IMPORTED and ALIAS targets compile nothing and are skipped. Comments do not count, and
# line breaks, spacing and the letter case of the command names do not matter. This check reads the text of
# CMakeLists.txt only. What it could not follow is refused: a bracket comment, and add_subdirectory, whose targets
# are declared in another file.
if(txt MATCHES "#\\[=*\\[")
  message(FATAL_ERROR "CMakeLists.txt has a bracket comment (#[[ ... ]]), which this check does not read; use line comments")
endif()
string(REGEX REPLACE "#[^\n]*" "" code "${txt}")           # without the comments
string(REGEX REPLACE "[ \t\r\n;]+" " " code "${code}")     # on one line, with single spaces
string(REGEX REPLACE " ?([()]) ?" "\\1" code "${code}")    # and none next to a parenthesis
foreach(command add_executable add_library add_subdirectory target_compile_options)   # in lower case
  set(any_case "")
  string(LENGTH "${command}" length)
  math(EXPR last "${length} - 1")
  foreach(i RANGE ${last})
    string(SUBSTRING "${command}" ${i} 1 letter)
    string(TOUPPER "${letter}" capital)
    string(APPEND any_case "[${letter}${capital}]")
  endforeach()
  string(REGEX REPLACE "${any_case}\\(" "${command}(" code "${code}")
endforeach()
if(code MATCHES "add_subdirectory\\(")
  message(FATAL_ERROR "CMakeLists.txt uses add_subdirectory: this check does not see the targets declared there, and has to be extended first")
endif()
# the name, and the keywords that follow it (STATIC, OBJECT, INTERFACE, IMPORTED, ALIAS, EXCLUDE_FROM_ALL, ...)
string(REGEX MATCHALL "add_(executable|library)\\([^ ()]+( [A-Z_0-9]+)*" declarations "${code}")
if(NOT declarations)
  message(FATAL_ERROR "no add_executable or add_library found in CMakeLists.txt: this check would pass on anything")
endif()
foreach(declaration IN LISTS declarations)
  if(declaration MATCHES " (INTERFACE|IMPORTED|ALIAS)( |$)")
    continue()
  endif()
  string(REGEX REPLACE "^add_[a-z]+\\(([^ ()]+).*$" "\\1" target "${declaration}")
  string(FIND "${code}" "target_compile_options(${target} PRIVATE \${QPUFORGE_COMPILE_OPTIONS})" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "${target} is not built with QPUFORGE_COMPILE_OPTIONS: CMakeLists.txt needs the line target_compile_options(${target} PRIVATE \${QPUFORGE_COMPILE_OPTIONS})")
  endif()
  message(STATUS "${target}: QPUFORGE_COMPILE_OPTIONS")
endforeach()
