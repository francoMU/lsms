

add_test(NAME li_extended_test
        COMMAND
        "$<TARGET_FILE:lsms> < i_lsms > stdout"
        WORKING_DIRECTORY
        ${CMAKE_BINARY_DIR}/tests)


add_custom_target(files)
add_custom_command(
        TARGET files
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/i_lsms"
        "${CMAKE_BINARY_DIR}/tests/i_lsms")

add_custom_command(
        TARGET files
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_CURRENT_SOURCE_DIR}/v_li.0"
        "${CMAKE_BINARY_DIR}/tests/v_li.0")

add_dependencies(files lsms)
