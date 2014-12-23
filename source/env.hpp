#ifndef QM_ENV_HPP
#define QM_ENV_HPP

// Encapsulates program environment

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

#if OS == OS_LINUX
#	include <GL/glx.h>
#	include <time.h>
#	include <X11/X.h>
#	include <X11/Xlib.h>
#elif OS == OS_WINDOWS
#	include <Windows.h>
#elif OS == OS_OSX // Using SDL on OSX because the native way is insane.
#	include <SDL.h>
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

#if OS == OS_LINUX
	Display* dpy;
	Window win;
	GLXContext ctx;
	timespec ts;
#elif OS == OS_WINDOWS
	HDC hDC;
	HWND hWnd;
	HGLRC hGlrc;
	DWORD ticks;
	static bool closeEvent;
#elif OS == OS_OSX
	SDL_Window* win;
	SDL_GLContext ctx;
	unsigned int ticks;
#endif // OS == OS_WINDOWS
};

Env envInit();
void envQuit(Env& env);
void envUpdate(Env& env);

typedef void (*voidFunc)();
voidFunc queryGlFunc(const char* name);

} // qm

#endif // QM_ENV_HPP
