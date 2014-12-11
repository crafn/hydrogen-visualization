// See unity.cpp for build instructions

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "env.hpp"
#include "gl.hpp"
#include "gui.hpp"

namespace qm {
struct Shader {
	GLuint vs;
	GLuint fs;
	GLuint prog;
	GLuint phaseLoc;
	GLuint colorLoc;
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
uniform vec3 u_color; \
varying vec3 v_pos; \
varying vec3 v_normal; \
/* x emission, y absorption */ \
vec2 sample(vec3 p) \
{ \
	p.z += 0.2; \
	float beam_e= \
		abs(cos(p.z*10.0 - sign(p.z)*u_phase)*0.5 + 1.0)* \
			0.001/(p.x*p.x + p.y*p.y + 0.001); \
	float disc_e= \
		0.01/((p.z*p.z + 0.01)*(p.x*p.x + p.y*p.y)*100.0 + 0.1); \
	float hole_a= pow(min(1.0, 0.1/(dot(p, p))), 10.0); \
	float a= clamp(hole_a, 0.0, 1.0); \
    return vec2((disc_e + beam_e)*(1 - a), a + disc_e*7.0)*20.0; \
} \
void main() \
{ \
	vec3 n= normalize(v_normal); \
	vec3 color= u_color; \
	float intensity= 0.0; \
	const int steps= 45; \
	const float dl= 2.0/steps; \
	for (int i= 0; i < steps; ++i) { \
		vec2 s= sample(v_pos + n*2.0*float(steps - i - 1)/steps); \
		intensity= max(0, intensity + (s.x - s.y*intensity)*dl); \
	} \
	gl_FragColor= vec4(color*intensity, 1.0); \
} \
\n";

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
			shd.colorLoc= glGetUniformLocation(shd.prog, "u_color");
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

void frame(const Env& env, const Shader& shd, const Vec2f cursor)
{
	static float setting_r= 0.3;
	static float setting_g= 0.8;
	static float setting_b= 1.0;


	glClear(GL_COLOR_BUFFER_BIT);

	{ // Volume
		static float phase;
		static Vec2f prev_smooth;

		// Smooth rotating
		Vec2f smooth= prev_smooth*0.5 + cursor*0.5;
		prev_smooth= smooth;

		/// @todo *dt
		phase += 0.3;

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(smooth.length()*200.0, smooth.y, -smooth.x, 0.0);

		glUseProgram(shd.prog);
		glUniform1f(shd.phaseLoc, phase);
		glUniform3f(shd.colorLoc, setting_r, setting_g, setting_b);

		glBegin(GL_QUADS);
			glVertex3f(-1, -1, 1);
			glVertex3f(1, -1, 1);
			glVertex3f(1,  1, 1);
			glVertex3f(-1, 1, 1);
		glEnd();
	}

	{ // Gui
		struct Slider {
			const char* title;
			float min;
			float max;
			float& value;
			int decimals;
			bool hover;

			static float height() { return 0.05; }
			static float width() { return 0.4; }
			static float top(std::size_t i) { return 1.0 - i*height(); }
			static float bottom(std::size_t i) { return top(i + 1); }
			float fraction() const
			{ return (value - min)/(max - min); }
		};

		static Slider sliders[]= {
			{ "r",         0.0,  1.0, setting_r,  2, false },
			{ "g",         0.0,  1.0, setting_g,  2, false },
			{ "b",         0.0,  1.0, setting_b,  2, false },
		};
		const std::size_t slider_count= sizeof(sliders)/sizeof(*sliders);

		// User interaction
		for (std::size_t i= 0; i < slider_count; ++i) {
			Slider& s= sliders[i];
			if (	cursor.x + 1.0 >= 0.0 && cursor.x + 1.0 < s.width() &&
					cursor.y >= s.bottom(i) && cursor.y < s.top(i)) {
				if (env.lmbDown)
					s.value= (1.0 + cursor.x)/s.width()*(s.max - s.min);
				s.hover= true;
			} else {
				s.hover= false;
			}
		}

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glUseProgram(0);

		glBegin(GL_QUADS);
			// Background
			glColor4f(0.0, 0.0, 0.0, 0.2);
			glVertex2f(-1.0, 1.0);
			glVertex2f(-1.0, 1.0 - slider_count*Slider::height());
			glVertex2f(Slider::width() - 1.0, 1.0 - slider_count*Slider::height());
			glVertex2f(Slider::width() - 1.0, 1.0);

			for (std::size_t i= 0; i < slider_count; ++i) {
				Slider& s= sliders[i];
				float top= s.top(i);
				float bottom= s.bottom(i);
				float width= s.width()*s.fraction();

				if (!s.hover)
					glColor4f(0.2, 0.2, 0.2, 0.5);
				else
					glColor4f(1.0, 1.0, 1.0, 0.5);

				glVertex2f(-1.0, top);
				glVertex2f(-1.0, bottom);
				glVertex2f(-1.0 + width, bottom);
				glVertex2f(-1.0 + width, top);
			}
		glEnd();
	}
} 

} // qm

int main()
{
	qm::Env env;
	qm::Shader shd;
	qm::init(env, shd);

	while (!env.quitRequested) {
		qm::envUpdate(env);
		glViewport(0, 0, env.winSize.x, env.winSize.y);
		qm::frame(	env,
					shd,
					env.cursorPos);
	}

	qm::quit(env, shd);
}

