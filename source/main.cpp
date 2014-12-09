#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <GL/gl.h>

#include <X11/X.h>
#include <X11/Xlib.h>
#include <GL/glx.h>
#include <unistd.h>


// Required GL 2.1 funcs
typedef GLuint (*GlCreateShader)(GLenum);
GlCreateShader glCreateShader;
typedef void (*GlShaderSource)(GLuint, GLsizei, const GLchar**, const GLint*);
GlShaderSource glShaderSource;
typedef void (*GlCompileShader)(GLuint);
GlCompileShader glCompileShader;
typedef GLuint (*GlCreateProgram)();
GlCreateProgram glCreateProgram;
typedef void (*GlAttachShader)(GLuint, GLuint);
GlAttachShader glAttachShader;
typedef void (*GlLinkProgram)(GLuint);
GlLinkProgram glLinkProgram;
typedef void (*GlUseProgram)(GLuint);
GlUseProgram glUseProgram;
typedef void (*GlGetShaderiv)(GLuint, GLenum, GLint*);
GlGetShaderiv glGetShaderiv;
typedef void (*GlGetProgramiv)(GLuint, GLenum, GLint*);
GlGetProgramiv glGetProgramiv;
typedef void (*GlGetShaderInfoLog)(GLuint, GLsizei, GLsizei*, GLchar*);
GlGetShaderInfoLog glGetShaderInfoLog;
typedef void (*GlGetProgramInfoLog)(GLuint, GLsizei, GLsizei*, GLchar*);
GlGetProgramInfoLog glGetProgramInfoLog;
typedef void (*GlDetachShader)(GLuint, GLuint);
GlDetachShader glDetachShader;
typedef void (*GlDeleteShader)(GLuint);
GlDeleteShader glDeleteShader;
typedef void (*GlDeleteProgram)(GLuint);
GlDeleteProgram glDeleteProgram;

namespace qm {

struct Env {
	Display*	dpy;
	Window      win;
	GLXContext  ctx;
};

struct Shader {
	GLuint vs;
	GLuint fs;
	GLuint prog;
};

const GLchar* g_vs= "\
#version 120\n\
varying vec3 v_pos; \
void main() \
{ \
    gl_FrontColor= gl_Color; \
    gl_TexCoord[0]= gl_MultiTexCoord0; \
    gl_Position= gl_Vertex; \
	v_pos= (gl_ModelViewProjectionMatrix*gl_Vertex).xyz; \
} \
\n";

const GLchar* g_fs= "\
#version 120\n\
varying vec3 v_pos; \
void main() \
{ \
	float length= dot(v_pos, v_pos); \
    gl_FragColor= vec4(sin(v_pos.x*10.0)*0.5 + 0.5, cos(v_pos.y*10.0 + v_pos.x*5.0)*0.5 + 0.5, v_pos.z, sin(v_pos.x*10.0) + 1.0 - length); \
} \
\n";

typedef void (*voidFunc)();
voidFunc queryGlFunc(const char* name)
{
	voidFunc f=	glXGetProcAddressARB((const GLubyte*)name);
	if (!f) {
		std::printf("Failed to query function: %s", name);
		std::abort();
	}
	return f;
}

void checkShaderStatus(GLuint shd)
{
	GLint status;
	glGetShaderiv(shd, GL_COMPILE_STATUS, &status);
	if (!status) {
		const GLsizei max_len= 512;
		GLchar log[max_len];
		glGetShaderInfoLog(shd, max_len, NULL, log);
		std::printf("Shader compilation failed: %s", log);
		std::abort();
	}
}

void checkProgramStatus(GLuint prog)
{
	GLint link_status;
	glGetProgramiv(prog, GL_LINK_STATUS, &link_status);
	if (!link_status) {
		const GLsizei size= 512;
		GLchar log[size];
		glGetProgramInfoLog(prog, size, NULL, log);
		std::printf("Program link failed: %s", log);
		std::abort();
	}
}

void init(Env& env, Shader& shd)
{
	{ // Create window
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
							0, 0, 600, 600, 0,
							vi->depth,
							InputOutput,
							vi->visual,
							CWColormap | CWEventMask,
							&swa);
		XMapWindow(env.dpy, env.win);
		XStoreName(env.dpy, env.win, "QM Test");

		env.ctx= glXCreateContext(env.dpy, vi, NULL, GL_TRUE);
		glXMakeCurrent(env.dpy, env.win, env.ctx);
	}

	{ // Query necessary GL functions
		glCreateShader= (GlCreateShader)queryGlFunc("glCreateShader");
		glShaderSource= (GlShaderSource)queryGlFunc("glShaderSource");
		glCompileShader= (GlCompileShader)queryGlFunc("glCompileShader");
		glCreateProgram= (GlCreateProgram)queryGlFunc("glCreateProgram");
		glAttachShader= (GlAttachShader)queryGlFunc("glAttachShader");
		glLinkProgram= (GlLinkProgram)queryGlFunc("glLinkProgram");
		glUseProgram= (GlUseProgram)queryGlFunc("glUseProgram");
		glGetShaderiv= (GlGetShaderiv)queryGlFunc("glGetShaderiv");
		glGetProgramiv= (GlGetProgramiv)queryGlFunc("glGetProgramiv");
		glGetShaderInfoLog= (GlGetShaderInfoLog)queryGlFunc("glGetShaderInfoLog");
		glGetProgramInfoLog= (GlGetProgramInfoLog)queryGlFunc("glGetProgramInfoLog");
		glDetachShader= (GlDetachShader)queryGlFunc("glDetachShader");
		glDeleteShader= (GlDeleteShader)queryGlFunc("glDeleteShader");
		glDeleteProgram= (GlDeleteProgram)queryGlFunc("glDeleteProgram");
	}

	{ // Create shaders
		{ // Vertex
			shd.vs= glCreateShader(GL_VERTEX_SHADER);
			glShaderSource(shd.vs, 1, &g_vs, NULL);
			glCompileShader(shd.vs);
			checkShaderStatus(shd.vs);
		}
		{ // Fragment
			shd.fs= glCreateShader(GL_FRAGMENT_SHADER);
			glShaderSource(shd.fs, 1, &g_fs, NULL);
			glCompileShader(shd.fs);
			checkShaderStatus(shd.fs);
		}
		{ // Shader program
			shd.prog= glCreateProgram();
			glAttachShader(shd.prog, shd.vs);
			glAttachShader(shd.prog, shd.fs);
			glLinkProgram(shd.prog);
			checkProgramStatus(shd.prog);
		}
	}

	{ // Setup initial GL state
		glClearColor(0.0, 0.0, 0.0, 1.0);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
}

void quit(const Env& env, const Shader& shd)
{
	{ // Close window
		glXMakeCurrent(env.dpy, None, NULL);
		glXDestroyContext(env.dpy, env.ctx);
		XDestroyWindow(env.dpy, env.win);
		XCloseDisplay(env.dpy);
	}

	{ // Destroy shaders
		glDetachShader(shd.prog, shd.vs);
		glDeleteShader(shd.vs);

		glDetachShader(shd.prog, shd.fs);
		glDeleteShader(shd.fs);

		glDeleteProgram(shd.prog);
	}
}

void draw(const Shader& shd)
{
	glUseProgram(shd.prog);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	static float rot;
	rot += 5;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(rot, 1.0, 0.1, 1.0);

	float step= 0.1;

	glBegin(GL_QUADS);
	for (int i= 0; i < 100; ++i) {
		glVertex3f(-.75, -.75, step*i);
		glVertex3f( .75, -.75, step*i);
		glVertex3f( .75,  .75, step*i);
		glVertex3f(-.75,  .75, step*i);
	}
	glEnd();
} 

bool loop(const Env& env, const Shader& shd)
{
	XWindowAttributes gwa;
	XGetWindowAttributes(env.dpy, env.win, &gwa);
	glViewport(0, 0, gwa.width, gwa.height);

	qm::draw(shd);

	usleep(1);
	glXSwapBuffers(env.dpy, env.win);

	while(XPending(env.dpy)) {
		XEvent xev;
        XNextEvent(env.dpy, &xev);
		if(xev.type == KeyPress)
			return false;
	}
	
	return true;
}

} // qm

int main()
{
	qm::Env env;
	qm::Shader shd;

	qm::init(env, shd);
	while(loop(env, shd));
	qm::quit(env, shd);
}

