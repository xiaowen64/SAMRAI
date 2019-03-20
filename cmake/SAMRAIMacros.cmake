macro (samrai_add_tests)

  set(singleValueArgs NAME EXECUTABLE PARALLEL)
  set(multiValueArgs INPUTS)
  set(counter 0)
  set(base_name ${arg_NAME})

  cmake_parse_arguments(arg
      "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN} )

  foreach (test_file ${arg_INPUTS})

    math(EXPR counter "${counter}+1")


    message(STATUS "Test: ${arg_NAME} with input ${test_file}")

    get_filename_component(short_test_file ${test_file} NAME)
    set(test_name "${base_name}_test_${short_test_file}")

    blt_add_test(NAME ${test_name}
      COMMAND ${arg_EXECUTABLE} ${test_file}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    if(${arg_PARALLEL})
      set(test_name "${base_name}_test_${short_test_file}_2")
      blt_add_test(NAME ${test_name}
        COMMAND ${arg_EXECUTABLE} ${test_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        NUM_MPI_TASKS 2)

      set(test_name "${base_name}_test_${short_test_file}_4")
      blt_add_test(NAME ${test_name}
        COMMAND ${arg_EXECUTABLE} ${test_file}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        NUM_MPI_TASKS 4)
    endif()

    set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  endforeach ()
endmacro ()
