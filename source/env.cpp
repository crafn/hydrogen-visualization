#include <cstdlib>
#include <cstdio>
#include "env.hpp"

#if PLATFORM == PLATFORM_LINUX
#	include <unistd.h>
#endif

namespace qm {

#if PLATFORM == PLATFORM_WINDOWS
bool Env::closeEvent;
#endif

Env envInit()
{
	const char* title= "QM Test";
	Vec2i reso(800, 800);

	Env env;
	env.cursorPos= Vec2f(0, 0);
	env.lmbDown= false;
	env.dt= 0.0;
	env.winSize= reso;
	env.quitRequested= false;

#if PLATFORM == PLATFORM_LINUX
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
	swa.event_mask= ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask;
	env.win=
		XCreateWindow(	env.dpy,
						root,
						0, 0, reso.x, reso.y, 0,
						vi->depth,
						InputOutput,
						vi->visual,
						CWColormap | CWEventMask,
						&swa);
	XMapWindow(env.dpy, env.win);
	XStoreName(env.dpy, env.win, title);
	
	env.ctx= glXCreateContext(env.dpy, vi, NULL, GL_TRUE);
	glXMakeCurrent(env.dpy, env.win, env.ctx);

	clock_gettime(CLOCK_MONOTONIC, &env.ts);

#elif PLATFORM == PLATFORM_WINDOWS
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

	WNDCLASS wc= WNDCLASS(); 
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
		0, 0, reso.x, reso.y, 0, 0, wc.hInstance, 0);

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
	
	env.ticks= GetTickCount();

#elif PLATFORM == PLATFORM_SDL
	SDL_Init(SDL_INIT_VIDEO);
	env.win= SDL_CreateWindow(
		title,
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		reso.x,
		reso.y,
		SDL_WINDOW_OPENGL);
	if (!env.win) {
		std::printf("SDL Error: %s\n", SDL_GetError());
		std::abort();
	}

	env.ctx= SDL_GL_CreateContext(env.win);
	env.ticks= SDL_GetTicks();
#endif

	return env;
}

void envQuit(Env& env)
{
#if PLATFORM == PLATFORM_LINUX
	glXMakeCurrent(env.dpy, None, NULL);
	glXDestroyContext(env.dpy, env.ctx);
	XDestroyWindow(env.dpy, env.win);
	XCloseDisplay(env.dpy);
#elif PLATFORM == PLATFORM_WINDOWS
	wglDeleteContext(env.hGlrc);
#elif PLATFORM == PLATFORM_SDL
	SDL_GL_DeleteContext(env.ctx);
	SDL_DestroyWindow(env.win);
	SDL_Quit();
#endif
}

void envUpdate(Env& env)
{
	Vec2f prev_cursor_pos= env.cursorPos;

#if PLATFORM == PLATFORM_LINUX
	usleep(1);
	glXSwapBuffers(env.dpy, env.win);

	while(XPending(env.dpy)) {
		XEvent xev;
        XNextEvent(env.dpy, &xev);
		if(xev.type == KeyPress) {
			int keys_ret;
			KeySym* keysym=
				XGetKeyboardMapping(env.dpy, xev.xkey.keycode, 1, &keys_ret);
			
			if (*keysym == XK_Escape)
				env.quitRequested= true;

			XFree(keysym);
		}

		if (xev.xbutton.type == ButtonPress)
			env.lmbDown= true;
		else if (xev.xbutton.type == ButtonRelease)
			env.lmbDown= false;
	}

	XWindowAttributes gwa;
	XGetWindowAttributes(env.dpy, env.win, &gwa);
	env.winSize.x= gwa.width;
	env.winSize.y= gwa.height;

	int root_x= 0, root_y= 0;
	Window w;
	int cursor_x, cursor_y;
	unsigned int mask;
	XQueryPointer(	env.dpy, env.win, &w,
					&w, &root_x, &root_y, &cursor_x, &cursor_y,
					&mask);

	env.cursorPos.x= 2.0*cursor_x/gwa.width - 1.0;
	env.cursorPos.y= 1.0 - 2.0*cursor_y/gwa.height;

	long old_us= env.ts.tv_nsec/1000 + env.ts.tv_sec*1000000;
	clock_gettime(CLOCK_MONOTONIC, &env.ts);
	long new_us= env.ts.tv_nsec/1000 + env.ts.tv_sec*1000000;
	env.dt= (new_us - old_us)/1000000.0;

#elif PLATFORM == PLATFORM_WINDOWS
	Sleep(1);
	SwapBuffers(env.hDC);

	MSG msg;
	while(PeekMessage(&msg, env.hWnd, 0, 0, PM_REMOVE) > 0) { 
		TranslateMessage(&msg); 
		DispatchMessage(&msg); 
	}
	env.lmbDown= GetKeyState(VK_LBUTTON) & 0x8000;

	if (Env::closeEvent)
		env.quitRequested= true;

	RECT rect;
	if(GetClientRect(env.hWnd, &rect)) {
		env.winSize.x= rect.right - rect.left;
		env.winSize.y= rect.bottom - rect.top;
	}
	POINT cursor;
	GetCursorPos(&cursor);
	ScreenToClient(env.hWnd, &cursor);
	env.cursorPos.x= 2.0*cursor.x/(rect.right - rect.left) - 1.0;
	env.cursorPos.y= 1.0 - 2.0*cursor.y/(rect.bottom - rect.top);

	DWORD old_ticks= env.ticks;
	env.ticks= GetTickCount();
	DWORD new_ticks= env.ticks;
	env.dt= (new_ticks - old_ticks)/1000.0;

#elif PLATFORM == PLATFORM_SDL
	SDL_Delay(1);
	SDL_GL_SwapWindow(env.win);

	SDL_Event e;
	while(SDL_PollEvent(&e)) {
		if (e.type == SDL_QUIT)
			env.quitRequested= true;	
	}

	int x, y;
	int state= SDL_GetMouseState(&x, &y);

	env.cursorPos.x= 2.0*x/env.winSize.x - 1.0;
	env.cursorPos.y= 1.0 - 2.0*y/env.winSize.y;

	bool was_down= env.lmbDown;
	env.lmbDown= state & SDL_BUTTON(SDL_BUTTON_LEFT);
	if (!was_down && env.lmbDown) {
		env.anchorPos= env.cursorPos;
	}

	unsigned int old_ticks= env.ticks;
	env.ticks= SDL_GetTicks();
	env.dt= (old_ticks - env.ticks)/1000.0;

#endif

	env.cursorDelta= env.cursorPos - prev_cursor_pos;
	if (!env.lmbDown)
		env.anchorPos= env.cursorPos;
}

typedef void (*voidFunc)();
voidFunc queryGlFunc(const char* name)
{
	voidFunc f= NULL;
#if PLATFORM == PLATFORM_LINUX
	f= glXGetProcAddressARB((const GLubyte*)name);
#elif PLATFORM == PLATFORM_WINDOWS
	f= (voidFunc)wglGetProcAddress(name);
#elif PLATFORM == PLATFORM_SDL
	f= (voidFunc)SDL_GL_GetProcAddress(name);
#endif
	if (!f) {
		std::printf("Failed to query gl function: %s\n", name);
		std::abort();
	}
	return f;
}

} // qm
