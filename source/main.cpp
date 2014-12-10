#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
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
typedef GLuint (*GlGetUniformLocation)(GLuint, const GLchar*);
GlGetUniformLocation glGetUniformLocation;
typedef void (*GlUniform1f)(GLuint, GLfloat);
GlUniform1f glUniform1f;

namespace qm {

struct Env {
	Display*	dpy;
	Window		win;
	GLXContext	ctx;
};

struct Shader {
	GLuint vs;
	GLuint fs;
	GLuint prog;
	GLuint phaseLoc;
};

const GLchar* g_vs= "\
#version 120\n\
varying vec3 v_pos; \
varying vec3 v_normal; \
void main() \
{ \
    gl_Position= vec4(gl_Vertex.xy, 0.0, 1.0); \
	v_pos= (gl_ModelViewMatrix*gl_Vertex).xyz; \
	v_normal= mat3(gl_ModelViewMatrix)*vec3(gl_Vertex.xy, -1.0); \
} \
\n";

const GLchar* g_fs= "\
#version 120\n\
uniform float u_phase; \
varying vec3 v_pos; \
varying vec3 v_normal; \
/* x emission, y translucency */ \
vec2 sample(vec3 p) \
{ \
	float beam_a= \
		abs(cos(p.x*10.0 - sign(p.x)*u_phase)*0.5 + 1.0)* \
			0.001/(p.z*p.z + p.y*p.y) + \
		0.01/((p.x*p.x + 0.01)*(p.z*p.z + p.y*p.y)*100.0 + 0.1); \
	float hole_a= pow(min(1.0, 0.1/(dot(p, p))), 10.0); \
	float lerp= clamp((1.0 - hole_a), 0.0, 1.0); \
    return vec2(beam_a*lerp, lerp); \
} \
void main() \
{ \
	vec3 n= normalize(v_normal); \
	vec3 color= vec3(0.3, 0.8, 1.0); \
	vec3 accum= vec3(0, 0, 0); \
	float transparency= 1.0; \
	const int steps= 45; \
	for (int i= 0; i < steps; ++i) { \
		vec2 s= sample(v_pos + n*2.0*float(i)/steps); \
		accum += color*s.x*transparency; \
		transparency *= s.y; \
	} \
	gl_FragColor= vec4(accum, 1.0); \
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
		glGetUniformLocation= (GlGetUniformLocation)queryGlFunc("glGetUniformLocation");
		glUniform1f= (GlUniform1f)queryGlFunc("glUniform1f");
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

			shd.phaseLoc= glGetUniformLocation(shd.prog, "u_phase");
		}
	}

	{ // Setup initial GL state
		glClearColor(0.0, 0.0, 0.0, 0.0);
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

void draw(const Shader& shd, float x, float y)
{
	static float phase;
	static float prev_x, prev_y;
	phase += 0.3;

	// Smooth rotating
	x= prev_x*0.5 + x*0.5;
	y= prev_y*0.5 + y*0.5;
	prev_x= x;
	prev_y= y;

	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(std::sqrt(x*x + y*y)*200.0, y, -x, 0.0);

	glUseProgram(shd.prog);
	glUniform1f(shd.phaseLoc, phase);

	glBegin(GL_QUADS);
		glVertex3f(-1, -1, 1);
		glVertex3f(1, -1, 1);
		glVertex3f(1,  1, 1);
		glVertex3f(-1, 1, 1);
	glEnd();
} 

bool loop(const Env& env, const Shader& shd)
{
	int root_x= 0, root_y= 0;
	Window w;
	int win_x, win_y;
	unsigned int mask_return;
	XQueryPointer(	env.dpy, env.win, &w,
					&w, &root_x, &root_y, &win_x, &win_y,
					&mask_return);

	XWindowAttributes gwa;
	XGetWindowAttributes(env.dpy, env.win, &gwa);
	glViewport(0, 0, gwa.width, gwa.height);
	
	qm::draw(	shd,
				2.0*win_x/gwa.width - 1.0,
				1.0 - 2.0*win_y/gwa.height);

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

