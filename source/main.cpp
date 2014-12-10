// See unity.cpp for build instructions

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "env.hpp"
#include "gl.hpp"

namespace qm {

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
			0.001/(p.z*p.z + p.y*p.y); \
	float disc_a= \
		0.01/((p.x*p.x + 0.01)*(p.z*p.z + p.y*p.y)*100.0 + 0.1); \
	float hole_a= pow(min(1.0, 0.1/(dot(p, p))), 10.0); \
	float lerp= clamp((1.0 - hole_a), 0.0, 1.0); \
    return vec2((disc_a + beam_a)*lerp, lerp); \
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
	env= envInit();
	queryGlFuncs();

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

void quit(Env& env, const Shader& shd)
{
	glDetachShader(shd.prog, shd.vs);
	glDeleteShader(shd.vs);

	glDetachShader(shd.prog, shd.fs);
	glDeleteShader(shd.fs);

	glDeleteProgram(shd.prog);
	
	envQuit(env);
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

} // qm

int main()
{
	qm::Env env;
	qm::Shader shd;
	qm::init(env, shd);

	while (!env.quitRequested) {
		envUpdate(env);
		glViewport(0, 0, env.winWidth, env.winHeight);
		qm::draw(	shd,
					env.cursorX,
					env.cursorY);
	}

	qm::quit(env, shd);
}

