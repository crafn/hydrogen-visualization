// See unity.cpp for build instructions

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "env.hpp"
#include "fontdata.hpp"
#include "gl.hpp"

namespace qm {

struct Program {
	struct VolumeShader {
		GLuint vs;
		GLuint fs;
		GLuint prog;
		GLuint phaseLoc;
		GLuint colorLoc;
	};
	struct VolumeFbo {
		GLuint fboId;
		GLuint texId;
		Vec2i reso;
		bool filtering;
	};
	struct Font {
		GLuint texId;
		Vec2f uv[256]; // Lower left corners of characters
		Vec2f charUvSize;
	};

	VolumeShader shader;
	VolumeFbo fbo;
	Font font;

	float time;
	float sampleCount;
	float resoMul;
	float filtering;
	float r, g, b;
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
uniform float u_phase; \
uniform vec3 u_color; \
varying vec3 v_pos; \
varying vec3 v_normal; \
float rand(vec2 co){ \
    return fract(sin(dot(co.xy, vec2(12.9898,78.233)))*43758.5453); \
} \
/* x emission, y absorption */ \
vec2 sample(vec3 p) \
{ \
	p.z += 0.2; \
	float beam_e= \
		abs(cos(p.z*10.0 - sign(p.z)*u_phase*10.0)*0.5 + 1.0)* \
			0.001/(p.x*p.x + p.y*p.y + 0.001); \
	float disc_e= \
		0.01/((p.z*p.z + 0.01)*(p.x*p.x + p.y*p.y)*100.0 + 0.1); \
	float hole_a= pow(min(1.0, 0.1/(dot(p, p))), 10.0); \
	float a= clamp(hole_a, 0.0, 1.0); \
    return vec2((disc_e + beam_e)*(1 - a), a + disc_e*7.0)*25.0; \
} \
void main() \
{ \
	vec3 n= normalize(v_normal); \
	vec3 color= u_color; \
	float intensity= 0.0; \
	const int steps= SAMPLE_COUNT; \
	const float dl= 2.0/steps; \
	for (int i= 0; i < steps; ++i) { \
		vec2 s= sample(v_pos + n*2.0*float(steps - i - 1)/steps); \
		intensity= max(0, intensity + (s.x - s.y*intensity)*dl); \
	} \
	intensity += rand(v_pos.xy + vec2(u_phase/10000.0, 0))*0.02; \
	gl_FragColor= vec4(color*intensity, 1.0); \
} \
\n";

Program::VolumeShader createShader(int sample_count)
{
	Program::VolumeShader shd;
	{ // Vertex
		shd.vs= glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(shd.vs, 1, &g_vs, NULL);
		glCompileShader(shd.vs);
		checkShaderStatus(shd.vs);
	}
	{ // Fragment
		const std::size_t buf_size= 512;
		char buf[buf_size];
		std::sprintf(buf,
			"#version 120\n"
			"#define SAMPLE_COUNT %i\n", sample_count);
		const GLchar* shader_src[]= { buf, g_fs };
		const GLsizei shader_src_count= sizeof(shader_src)/sizeof(*shader_src);

		shd.fs= glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(shd.fs, shader_src_count, shader_src, NULL);
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
	return shd;
}

void destroyShader(Program::VolumeShader& shd)
{
	glDetachShader(shd.prog, shd.vs);
	glDeleteShader(shd.vs);

	glDetachShader(shd.prog, shd.fs);
	glDeleteShader(shd.fs);

	glDeleteProgram(shd.prog);
}

Program::VolumeFbo createFbo(Vec2i reso, bool filtering)
{
	GLenum filter= filtering ? GL_LINEAR : GL_NEAREST;

	Program::VolumeFbo fbo;
	fbo.reso= reso;
	fbo.filtering= filtering;
	glGenTextures(1, &fbo.texId);
	glBindTexture(GL_TEXTURE_2D, fbo.texId);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexImage2D(	GL_TEXTURE_2D, 0, GL_RGB,
					fbo.reso.x, fbo.reso.y,
					0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

	glGenFramebuffers(1, &fbo.fboId);
	glGenFramebuffers(1, &fbo.fboId);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo.fboId);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo.texId, 0);
	return fbo;
}

void destroyFbo(Program::VolumeFbo& fbo)
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glDeleteFramebuffers(1, &fbo.fboId);
	glDeleteTextures(1, &fbo.texId);
}

void init(Env& env, Program& prog)
{
	env= envInit();
	queryGlFuncs();

	{ // Font
		Program::Font& font= prog.font;
		Vec2f char_size= cast<Vec2f>(g_font.charSize);
		Vec2f font_size= cast<Vec2f>(g_font.size);
		font.charUvSize= char_size/font_size; 

		Vec2i cursor(0, g_font.size.y - g_font.charSize.y); // Lower-left origin
		for (std::size_t i= 0; i < std::strlen(g_font.chars); ++i) {
			unsigned char ch= g_font.chars[i];
			font.uv[ch]= cast<Vec2f>(cursor)/font_size;

			cursor.x += g_font.charSize.x;
			if (cursor.x + g_font.charSize.x > g_font.size.x) {
				cursor.x= 0;
				cursor.y -= g_font.charSize.y;
			}
		}

		glGenTextures(1, &font.texId);
		glBindTexture(GL_TEXTURE_2D, font.texId);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(	GL_TEXTURE_2D, 0, GL_ALPHA8,
						g_font.size.x, g_font.size.y,
						0, GL_ALPHA, GL_UNSIGNED_BYTE,
						g_font.data);
	}

	{ // Program state
		prog.time= 0.0;
		prog.sampleCount= 40;
		prog.resoMul= 0.5;
		prog.filtering= 0.0;
		prog.r= 0.4;
		prog.g= 0.8;
		prog.b= 1.0;

		prog.shader= createShader(prog.sampleCount);
		prog.fbo= createFbo(env.winSize*prog.resoMul, prog.filtering > 0.5);
	}

	{ // Setup initial GL state
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
}

void quit(Env& env, Program& prog)
{
	destroyFbo(prog.fbo);
	destroyShader(prog.shader);

	{ // Font
		Program::Font& font= prog.font;
		glDeleteTextures(1, &font.texId);
	}

	envQuit(env);
}

void frame(const Env& env, Program& prog)
{
	struct Slider {
		const char* title;
		float min;
		float max;
		float& value;
		int decimals;
		bool recompile;

		static float height() { return 0.05; }
		static float width() { return 0.65; }
		static float top(std::size_t i) { return 1.0 - i*height(); }
		static float bottom(std::size_t i) { return top(i + 1); }
		bool pointInside(std::size_t i, Vec2f p) const
		{
			return	p.x  >= -1.0	&& p.x < width() - 1.0 &&
					p.y > bottom(i)	&& p.y < top(i);
		}
		float fraction() const
		{
			float f= (value - min)/(max - min);
			return f < 1 ? f : 1;
		}
		float coordToValue(float x) const
		{
			float v= clamp((1.0 + x)/width()*(max - min) + min, min, max);
			return round(v, decimals);
		}
	};

	static Slider sliders[]= {
		{ "Time",		0.0,	5.0,	prog.time,			6, false },
		{ "Samples",	5,		150,	prog.sampleCount,	0, true },
		{ "Resolution",	0.01,	1.0,	prog.resoMul,		6, false },
		{ "Filtering",	0,		1,		prog.filtering,		0, false },
		{ "R",			0.0,	2.0,	prog.r,				6, false },
		{ "G",			0.0,	2.0,	prog.g,				6, false },
		{ "B",			0.0,	2.0,	prog.b,				6, false },
	};
	const std::size_t slider_count= sizeof(sliders)/sizeof(*sliders);
	bool slider_hover[slider_count]= {};

	static Vec2f prev_delta;
	static Vec2f rot;

	{ // User interaction
		bool slider_activity= false;
		for (std::size_t i= 0; i < slider_count; ++i) {
			Slider& s= sliders[i];
			if (s.pointInside(i, env.anchorPos)) {
				bool value_changed= false;
				if (env.lmbDown) {
					float new_value= s.coordToValue(env.cursorPos.x);
					value_changed= new_value != s.value;
					s.value= new_value;
				}
				slider_hover[i]= true;
				slider_activity= true;

				if (value_changed && s.recompile) {
					destroyShader(prog.shader);
					prog.shader= createShader(prog.sampleCount);
				}
			} else {
				slider_hover[i]= false;
			}
		}

		if (env.lmbDown && !slider_activity) {
			Vec2f smooth_delta= prev_delta*0.5 + (env.cursorDelta)*0.5;
			prev_delta= smooth_delta;
			
			rot += smooth_delta*2;
			rot.y= clamp(rot.y, -tau/4, tau/4);
		}
	}

	{ // Adjust FBO
		Vec2i volume_reso= cast<Vec2i>(cast<Vec2f>(env.winSize)*prog.resoMul);
		bool volume_filtering= prog.filtering > 0.5;
		if (volume_reso != prog.fbo.reso || volume_filtering != prog.fbo.filtering) {
			destroyFbo(prog.fbo);
			prog.fbo= createFbo(volume_reso, volume_filtering);
		}
	}

	/// @todo dt
	prog.time += 1.0/30.0;

	glClear(GL_COLOR_BUFFER_BIT);

	/// @todo Fixed pipeline starts to feel a bit clumsy
	{ // Draw volume	
		// Draw to fbo
		glBindFramebuffer(GL_FRAMEBUFFER, prog.fbo.fboId);
		glViewport(0, 0, prog.fbo.reso.x, prog.fbo.reso.y);
		glDisable(GL_TEXTURE_2D);

		Program::VolumeShader& shd= prog.shader;

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		// Trackball-style rotation
		glRotatef(rot.y*radToDeg, cos(rot.x), 0.0, sin(rot.x));
		glRotatef(rot.x*radToDeg, 0.0, -1.0, 0.0);

		glUseProgram(shd.prog);
		glUniform1f(shd.phaseLoc, prog.time);
		glUniform3f(shd.colorLoc, prog.r, prog.g, prog.b);

		glBegin(GL_QUADS);
			glVertex3f(-1, -1, 1);
			glVertex3f(1, -1, 1);
			glVertex3f(1,  1, 1);
			glVertex3f(-1, 1, 1);
		glEnd();

		// Draw scaled fbo texture
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, env.winSize.x, env.winSize.y);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, prog.fbo.texId);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glMatrixMode(GL_TEXTURE);
		glLoadIdentity();
		glUseProgram(0);
		glBegin(GL_QUADS);
			glTexCoord2f(0, 0); glVertex2f(-1, -1);
			glTexCoord2f(1, 0); glVertex2f(1, -1);
			glTexCoord2f(1, 1); glVertex2f(1, 1);
			glTexCoord2f(0, 1); glVertex2f(-1, 1);
		glEnd();
	}

	{ // Draw gui
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, env.winSize.x, env.winSize.y);
		glMatrixMode(GL_TEXTURE);
		glLoadIdentity();
		glScalef(1.0, -1.0, 1.0); // Textures are upside-down

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glUseProgram(0);
		glDisable(GL_TEXTURE_2D);

		glBegin(GL_QUADS);
			// Background
			glColor4f(0.1, 0.1, 0.1, 0.3);
			glVertex2f(-1.0, 1.0);
			glVertex2f(-1.0, 1.0 - slider_count*Slider::height());
			glVertex2f(Slider::width() - 1.0, 1.0 - slider_count*Slider::height());
			glVertex2f(Slider::width() - 1.0, 1.0);

			// Sliders
			for (std::size_t i= 0; i < slider_count; ++i) {
				Slider& s= sliders[i];
				float top= s.top(i);
				float bottom= s.bottom(i);
				float width= s.width()*s.fraction();

				if (!slider_hover[i])
					glColor4f(0.3, 0.3, 0.3, 0.6);
				else
					glColor4f(1.0, 1.0, 1.0, 0.8);

				glVertex2f(-1.0, top);
				glVertex2f(-1.0, bottom);
				glVertex2f(-1.0 + width, bottom);
				glVertex2f(-1.0 + width, top);
			}
		glEnd();

		// Slider texts
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, prog.font.texId);
		glBegin(GL_QUADS);
		glColor4f(0.8, 0.8, 0.8, 1.0);
		for (std::size_t s_i= 0; s_i < slider_count; ++s_i) {
			Slider& s= sliders[s_i];
			Vec2f pos= fitToGrid(Vec2f(-0.98, s.bottom(s_i)), env.winSize);

			for (std::size_t c_i= 0; c_i < std::strlen(s.title); ++c_i) {
				unsigned char ch= s.title[c_i];
				Vec2f ch_size= cast<Vec2f>(g_font.charSize)/cast<Vec2f>(env.winSize)*2.0;
				Vec2f ch_pos(pos.x + ch_size.x*c_i, pos.y);

				Vec2f ch_ur= ch_pos + ch_size;
				Vec2f ll= prog.font.uv[ch];
				Vec2f ur= ll + prog.font.charUvSize;
				glTexCoord2f(ll.x, ll.y); glVertex2f(ch_pos.x, ch_pos.y);
				glTexCoord2f(ur.x, ll.y); glVertex2f(ch_ur.x,  ch_pos.y);
				glTexCoord2f(ur.x, ur.y); glVertex2f(ch_ur.x,  ch_ur.y);
				glTexCoord2f(ll.x, ur.y); glVertex2f(ch_pos.x, ch_ur.y);
			}
		}
		glEnd();
	}

#ifndef NDEBUG
	checkGlErrors("frame end");
#endif
} 

} // qm

int main()
{
	qm::Env env;
	qm::Program prog;
	qm::init(env, prog);

	while (!env.quitRequested) {
		qm::envUpdate(env);
		qm::frame(env, prog);
	}

	qm::quit(env, prog);
}

