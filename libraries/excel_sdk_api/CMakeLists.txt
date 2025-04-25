project(excel_sdk_api)

set(includes includes/)
set(sources src/xlframework.c
            src/xlMemoryManager.cpp
            src/xlMemoryPool.cpp
            src/xlcall.cpp
)

add_library(${PROJECT_NAME} ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes} ${common_includes_dir})