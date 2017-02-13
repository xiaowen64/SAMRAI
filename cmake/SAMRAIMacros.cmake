macro (samrai_add_tests)

  set(singleValueArgs NAME EXECUTABLE PARALLEL)
  set(multiValueArgs INPUTS)
  set(counter 0)

  cmake_parse_arguments(arg
      "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

  foreach (test_file ${arg_INPUTS})

    math(EXPR counter "${counter}+1")

    message(STATUS "Test: ${arg_NAME} with input ${test_file}")

    get_filename_component(short_test_file ${test_file} NAME)
    set(test_name "${arg_NAME}_test_${short_test_file}")

    add_test(NAME ${test_name}
      COMMAND $<TARGET_FILE:${arg_EXECUTABLE}> ${test_file})

    set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  endforeach ()
endmacro ()
