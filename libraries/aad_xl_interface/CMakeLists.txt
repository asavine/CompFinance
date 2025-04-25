project(aad_xl_interface)

set(includes includes/ 
            ${aad_dir}/ 
            ${aad_dir}/includes 
            ${excel_sdk_api_dir}/ 
            ${excel_sdk_api_dir}/includes
            ${common_dir}/
)

set(sources src/xlExport.cpp)

add_compile_definitions(${PROJECT_NAME}_EXPORTS
                        _WINDLL
)

add_library(${PROJECT_NAME} SHARED ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes} ${common_includes_dir})
target_link_libraries(${PROJECT_NAME} PUBLIC aad excel_sdk_api)
set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".xll")