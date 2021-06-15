

add_library(mjson)

target_sources(mjson
        PUBLIC
        json.c
        json.h
        json_helper.c
        json_helper.h
        rstring.c
        rstring.h
        )

target_include_directories(
        mjson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)