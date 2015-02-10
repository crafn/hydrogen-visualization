//
// Building
//
// On Linux
// GCC: g++ -O2 source/unity.cpp -lGL -lX11 -o qm
// Clang: clang++ -O2 source/unity.cpp -lGL -lX11 -o qm
//
// On Windows
// MinGW: g++ -O2 source/unity.cpp -lOpenGL32 -lGdi32 -o qm.exe
// MSVC: TODO
//
// On Mac OS X
// Download and install SDL2 from https://www.libsdl.org/download-2.0.php
// Clang: clang -O2 source/unity.cpp -framework OpenGL -framework SDL2 -lstdc++ -o qm
//
// Native libraries are used on Windows and Linux, but it's possible to build using SDL also on
// these platforms if you have it installed, by adding -DPLATFORM=3 and -lSDL2 to the command

#include "env.cpp"
#include "fontdata.cpp"
#include "main.cpp"
