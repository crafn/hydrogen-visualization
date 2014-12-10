#include <cstdlib>
#include <cstdio>
#include "env.hpp"

#if OS == OS_LINUX
#	include <GL/glx.h>
#	include <unistd.h>
#endif

namespace qm {

#if OS == OS_WINDOWS
bool Env::closeEvent;
#endif

Env envInit()
{
	Env env;
	env.cursorX= 0.0;
	env.cursorY= 0.0;
	env.winWidth= 1;
	env.winHeight= 1;
	env.quitRequested= false;

	const char* title= "QM Test";
	int resox= 600;
	int resoy= 600;

#if OS == OS_LINUX
	env.dpy= XOpenDisplay(NULL);
	if(env.dpy == NULL)
		std::abort();

	Window root= DefaultRootWindow(env.dpy);
	GLint att[]= { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
	XVisualInfo* vi= glXChooseVisual(env.dpy, 0, att);

	if(vi == NULL)
		std::abort();

	Colormap cmap;
	cmap= XCreateColormap(env.dpy, root, vi->visual, AllocNone);
	XSetWindowAttributes swa;
	swa.colormap= cmap;
	swa.event_mask= ExposureMask | KeyPressMask;
	env.win=
		XCreateWindow(	env.dpy,
						root,
						0, 0, resox, resoy, 0,
						vi->depth,
						InputOutput,
						vi->visual,
						CWColormap | CWEventMask,
						&swa);
	XMapWindow(env.dpy, env.win);
	XStoreName(env.dpy, env.win, title);

	env.ctx= glXCreateContext(env.dpy, vi, NULL, GL_TRUE);
	glXMakeCurrent(env.dpy, env.win, env.ctx);

	XWindowAttributes gwa;
	XGetWindowAttributes(env.dpy, env.win, &gwa);
	env.winWidth= gwa.width;
	env.winHeight= gwa.height;

#elif OS == OS_WINDOWS
	struct WndProc {
		static LRESULT CALLBACK call(
			HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
		{
			switch (message) {
				case WM_DESTROY:
					PostQuitMessage(0);
				break;
				case WM_CLOSE:
					Env::closeEvent= true;
				break;
				case WM_MOUSEMOVE:
					SetCursor(LoadCursor(NULL, IDC_ARROW));
				break;
				default:
					return DefWindowProc(hWnd, message, wParam, lParam);
			}
			return 0;
		}
	};

	MSG msg= {0};
	WNDCLASS wc= {0}; 
	wc.lpfnWndProc= WndProc::call;
	wc.hInstance= GetModuleHandle(0);
	wc.hbrBackground= (HBRUSH)(COLOR_BACKGROUND);
	wc.lpszClassName= title;
	wc.style = CS_OWNDC;
	if( !RegisterClass(&wc) ) {
			std::printf("RegisterClass failed\n");
			std::abort();
	}
	env.hWnd= CreateWindow(
		wc.lpszClassName,
		title,
		WS_OVERLAPPEDWINDOW|WS_VISIBLE,
		0, 0, resox, resoy, 0, 0, wc.hInstance, 0);

	// Create OpenGL context
	PIXELFORMATDESCRIPTOR pfd= {
		sizeof(PIXELFORMATDESCRIPTOR), 1,
		PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,
		PFD_TYPE_RGBA,
		32, // Framebuffer
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		24, // Depth buffer
		8, // Stencil buffer
		0, PFD_MAIN_PLANE, 0, 0, 0, 0
	};

	env.hDC= GetDC(env.hWnd);
	int choose= ChoosePixelFormat(env.hDC, &pfd);
	SetPixelFormat(env.hDC, choose, &pfd);
	env.hGlrc= wglCreateContext(env.hDC);
	wglMakeCurrent(env.hDC, env.hGlrc);
#endif

	return env;
}

void envQuit(Env& env)
{
#if OS == OS_LINUX
	glXMakeCurrent(env.dpy, None, NULL);
	glXDestroyContext(env.dpy, env.ctx);
	XDestroyWindow(env.dpy, env.win);
	XCloseDisplay(env.dpy);
#elif OS == OS_WINDOWS
	wglDeleteContext(env.hGlrc);
#endif
}

void envUpdate(Env& env)
{
#if OS == OS_LINUX
	usleep(1);
	glXSwapBuffers(env.dpy, env.win);

	while(XPending(env.dpy)) {
		XEvent xev;
        XNextEvent(env.dpy, &xev);
		if(xev.type == KeyPress)
			env.quitRequested= true;
	}

	XWindowAttributes gwa;
	XGetWindowAttributes(env.dpy, env.win, &gwa);
	env.winWidth= gwa.width;
	env.winHeight= gwa.height;

	int root_x= 0, root_y= 0;
	Window w;
	int cursor_x, cursor_y;
	unsigned int mask;
	XQueryPointer(	env.dpy, env.win, &w,
					&w, &root_x, &root_y, &cursor_x, &cursor_y,
					&mask);

	env.cursorX= 2.0*cursor_x/gwa.width - 1.0;
	env.cursorY= 1.0 - 2.0*cursor_y/gwa.height;

#elif OS == OS_WINDOWS
	Sleep(1);
	SwapBuffers(env.hDC);
	
	MSG msg;
	BOOL bRet;
	while(PeekMessage(&msg, env.hWnd, 0, 0, PM_REMOVE) > 0) { 
		TranslateMessage(&msg); 
		DispatchMessage(&msg); 
	}

	if (Env::closeEvent)
		env.quitRequested= true;

	RECT rect;
	if(GetWindowRect(env.hWnd, &rect)) {
		env.winWidth= rect.right - rect.left;
		env.winHeight= rect.bottom - rect.top;
	}
	POINT cursor;
	GetCursorPos(&cursor);
	ScreenToClient(env.hWnd, &cursor);
	env.cursorX= 2.0*cursor.x/(rect.right - rect.left) - 1.0;
	env.cursorY= 1.0 - 2.0*cursor.y/(rect.bottom - rect.top);

#endif
}

typedef void (*voidFunc)();
voidFunc queryGlFunc(const char* name)
{
#if OS == OS_LINUX
	voidFunc f=	glXGetProcAddressARB((const GLubyte*)name);
#elif OS == OS_WINDOWS
	voidFunc f=	(voidFunc)wglGetProcAddress(name);
#endif
	if (!f) {
		std::printf("Failed to query gl function: %s\n", name);
		std::abort();
	}
	return f;
}

} // qm