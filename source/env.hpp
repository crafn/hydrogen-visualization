#ifndef QM_ENV_HPP
#define QM_ENV_HPP

// Encapsulates program environment

#define PLATFORM_WINDOWS 1
#define PLATFORM_LINUX 2
#define PLATFORM_SDL 3

#define OS_WINDOWS 1
#define OS_LINUX 2
#define OS_OSX 3

#ifdef __linux
#	define OS OS_LINUX
#elif _WIN32
#	define OS OS_WINDOWS
#elif __APPLE__
#	define OS OS_OSX
#endif

#if !defined(PLATFORM)
#	if OS == OS_LINUX
#		define PLATFORM PLATFORM_LINUX
#	elif OS == OS_WINDOWS
#		define PLATFORM PLATFORM_WINDOWS
#	elif OS == OS_OSX
#		define PLATFORM PLATFORM_SDL // Using SDL on OSX because the native way is insane.
#	endif
#endif

#if PLATFORM == PLATFORM_LINUX
#	include <GL/glx.h>
#	include <time.h>
#	include <X11/X.h>
#	include <X11/Xlib.h>
#	ifdef Complex
#		undef Complex
#	endif
#elif PLATFORM == PLATFORM_WINDOWS
#	include <Windows.h>
#elif PLATFORM == PLATFORM_SDL
#	include <SDL2/SDL.h>
#endif

#include "util.hpp"

namespace qm {

struct Env {
	Vec2f cursorPos; // Cursor position in OpenGL coordinates
	Vec2f anchorPos; // Mouse dragging start position
	Vec2f cursorDelta;
	bool lmbDown;
	float dt;
	Vec2i winSize; // Window content size in pixels
	bool quitRequested;

#if PLATFORM == PLATFORM_LINUX
	Display* dpy;
	Window win;
	GLXContext ctx;
	timespec ts;
#elif PLATFORM == PLATFORM_WINDOWS
	HDC hDC;
	HWND hWnd;
	HGLRC hGlrc;
	DWORD ticks;
	static bool closeMessage;
	static bool lbuttondownMessage;
#elif PLATFORM == PLATFORM_SDL
	SDL_Window* win;
	SDL_GLContext ctx;
	unsigned int ticks;
#endif
};

Env envInit();
void envQuit(Env& env);
void envUpdate(Env& env);

typedef void (*voidFunc)();
voidFunc queryGlFunc(const char* name);

} // qm

#endif // QM_ENV_HPP
