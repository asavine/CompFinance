project(aad)

set(includes includes/)

set(sources src/AAD.cpp
            src/mcBase.cpp
            src/sobol.cpp
            src/ThreadPool.cpp
)

# assuming at the moment that we are only msvc

add_library(${PROJECT_NAME} ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes} ${common_includes_dir})
