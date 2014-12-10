#ifndef QM_ENV_HPP
#define QM_ENV_HPP

// Encapsulates program environment

#define OS_WINDOWS 1
#define OS_LINUX 2
#ifdef __linux
#	define OS OS_LINUX
#elif _WIN32
#	define OS OS_WINDOWS
#endif

#if OS == OS_LINUX
#	include <X11/X.h>
#	include <X11/Xlib.h>
#elif OS == OS_WINDOWS
#	include <Windows.h>
#endif

namespace qm {

struct Env {
	// Cursor position in OpenGL coordinates
	float cursorX;
	float cursorY;
	
	// Window size in pixels
	int winWidth;
	int winHeight;
	
	bool quitRequested;

#if OS == OS_LINUX
	Display*	dpy;
	Window		win;
	GLXContext	ctx;
#elif OS == OS_WINDOWS
	HDC hDC;
	HWND hWnd;
	HGLRC hGlrc;
	static bool closeEvent;
#endif // OS == OS_WINDOWS
};

Env envInit();
void envQuit(Env& env);
void envUpdate(Env& env);

typedef void (*voidFunc)();
voidFunc queryGlFunc(const char* name);

} // qm

#endif // QM_ENV_HPP