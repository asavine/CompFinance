# Some Notes

Some useful commands, assuming we are in the root of the project:

For Windows:

1. mkdir build
2. cd build && cmake .. -G "Visual Studio 17 2022" -A x64
3. cmake --build . --config Release
4. cmake --build . --config Debug
5. cmake --install . -> this needs a pre setup if you are not admin, you could hit some errors.
