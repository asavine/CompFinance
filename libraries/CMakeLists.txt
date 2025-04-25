project(libraries)

set(common_includes_dir ${CMAKE_CURRENT_SOURCE_DIR}/common/)

add_subdirectory(aad)
set(aad_dir ${CMAKE_CURRENT_SOURCE_DIR}/aad)

if(WIN32)
    message(STATUS "Windows detected - xll interfaces")
    set(excel_sdk_api_dir ${CMAKE_CURRENT_SOURCE_DIR}/excel_sdk_api)

    add_subdirectory(excel_sdk_api)
    add_subdirectory(aad_xl_interface)
endif()
