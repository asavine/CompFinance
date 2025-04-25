@ECHO OFF

REM Defaults:
set toolset=v143
set config=Release
set my_cwd=%cd%
set clean=
set rebuild=N

echo "Current path %my_cwd%"

set /p toolset="Enter the toolset version for msvc [%toolset%] - only support for v141 - v142 and v143 with std c++20 >> "
set /p config="Enter the config to be built [%config%] - only support for Release or Debug >> "
set /p rebuild="Rebuild (ie cleaning before build) (Y/N)? [%rebuild%] >> "

if /I %rebuild% NEQ "N" (
	set clean="--clean-first"
)

set build_folder=build_%toolset%

if exist %build_folder% (
	cd %build_folder%
	cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_GENERATOR_TOOLSET=%toolset%
	cmake --build . --config %config% %clean%
) else (
	mkdir %build_folder%
	cd %build_folder%
	cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_GENERATOR_TOOLSET=%toolset%
	cmake --build . --config %config% %clean%
)

cd %my_cwd%
