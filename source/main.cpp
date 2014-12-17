// See unity.cpp for build instructions

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "env.hpp"
#include "fontdata.hpp"
#include "gl.hpp"

namespace qm {

struct Program {
	struct VolumeShader {
		GLuint vs, fs, prog;
		GLint timeLoc;
		GLint phaseLoc;
		GLint colorLoc;
		GLint transformLoc;
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
		Vec2f whiteTexelUv;
	};
	struct GuiShader {
		GLuint vs, fs, prog;
		GLint texLoc;
		GLint colorLoc;
	};
	struct QuadVbo {
		GLuint vboId;
	};

	VolumeShader shader;
	VolumeFbo fbo;
	Font font;
	GuiShader guiShader;
	QuadVbo vbo;
	float time;

	// Slider settings
	float phase;
	float sampleCount;
	float resoMul;
	float filtering;
	float r, g, b;
	float scale;
	float n, l, m;
};

Program::VolumeShader createVolumeShader(int sample_count, float scale, int n, int l, int m)
{
	/// @todo Remove
	testMath();

	if (l > n - 1)
		l= n - 1;
	if (m > l)
		m= l;

	double bohr_radius= std::exp(scale*0.5) - 0.99999;
	const int poly_term_count= 30;
	// Formula for hydrogen wave function with parameters r, theta, and phi
	StackString<1024> hydrogen_str;
	{
		// Form a hydrogen wave function |nlm> in four parts
		// psi_nlm(r, theta, phi) = C*E*L*Y, where
		//   C = normalization factor sqrt[(2/(n*a_0))^3*(n - l - 1)!/(2n(n + l)!)]
		//   E = e^(-rho/l)*rho^l
		//   L = Generalized Laguerre Polynomial L(n - l - 1, 2l + 1, rho)
		//   Y = Spherical harmonic function Y(l, m, theta, phi)
		//   rho = 2r/(n*a_0)
		const std::size_t rho_str_size= 16;
		char rho_str[rho_str_size];
		std::snprintf(rho_str, rho_str_size, "%e*r", 2.0/bohr_radius/n);

		// C
		hydrogen_str.append("%e",
			std::sqrt(	std::pow(2.0/(n*bohr_radius), 3)*
						fact(n - l - 1)/(2*n*fact(n + l))));
		hydrogen_str.append("*");

		// E
		hydrogen_str.append("exp(-%s/2.0)*pow(%s, %i.0)", rho_str, rho_str, l);
		hydrogen_str.append("*");

		{ // L
			hydrogen_str.append("(");
			double laguerre_coeff[poly_term_count]= {};
			assert(n - l < poly_term_count);

			laguerre(laguerre_coeff, n - l - 1, 2*l + 1);
			for (int i= 0; i < poly_term_count; ++i) {
				if (laguerre_coeff[i] == 0.0)
					continue;

				hydrogen_str.append("+%e*pow(%s, %i.0)", laguerre_coeff[i], rho_str, i);
			}
			hydrogen_str.append(")");
		}
		hydrogen_str.append("*");

		{ // Y
			// Dropping (-1)^m
			double normalization= 
				std::sqrt( (2*l + 1.0)/(4*pi)*fact(l - m)/fact(l + m) );
			hydrogen_str.append(
				"%e*pow(sin_theta, %i.0)*(",
				normalization,
				m); // Dropping e^(imphi)
			double spherical_coeff[poly_term_count]= {};
			assert(l - 1 < poly_term_count);
			sphericalHarmonics(spherical_coeff, l, m);
			for (int i= 0; i < poly_term_count; ++i) {
				if (spherical_coeff[i] == 0.0)
					continue;
				hydrogen_str.append("+%e*pow(cos_theta, %i)", spherical_coeff[i], i);
			}
			hydrogen_str.append(")");
		}

		std::printf("Wave function:\n%s\n", hydrogen_str.str);
	}

	const GLchar* vs_src=
		"#version 120\n"
		"attribute vec2 a_pos;"
		"attribute vec2 a_uv;"
		"uniform mat4 u_transform;"
		"varying vec3 v_pos;"
		"varying vec3 v_normal;"
		"void main()"
		"{"
		"	gl_Position= vec4(a_pos, 0.0, 1.0);"
		"	v_pos= (u_transform*vec4(a_pos, 0.0, 1.0)).xyz;"
		"	v_normal= mat3(u_transform)*vec3(a_pos, -1.0);"
		"}"
		"\n";
	const GLchar* fs_template_src=
		"uniform float u_phase;"
		"uniform float u_time;"
		"uniform vec3 u_color;"
		"varying vec3 v_pos;"
		"varying vec3 v_normal;"
		"float rand(vec2 co){"
		"    return fract(sin(dot(co.xy, vec2(12.9898,78.233)))*43758.5453);"
		"}"
		"vec2 sample(vec3 sphe_p)" // x emission, y absorption
		"{"
		"  float cos_theta= cos(sphe_p.y); /* @todo don't calculate */"
		"  float sin_theta= sin(sphe_p.y); /* @todo don't calculate */"
		"  float value= WAVEFUNC(sphe_p.x, sphe_p.y, sphe_p.z, cos_theta, sin_theta);"
		"  return vec2(value*value, 0.0);"
		"}"
		"void main()"
		"{"
		"	vec3 n= normalize(v_normal);"
		"	vec3 color= u_color;"
		"	float intensity= 0.0;"
		"	const float dl= 2.0/SAMPLE_COUNT;"
		"	for (int i= 0; i < SAMPLE_COUNT; ++i) {"
		"		vec3 cart_p= v_pos + n*2.0*float(SAMPLE_COUNT - i - 1)/SAMPLE_COUNT;"
		"		float dist= sqrt(dot(cart_p, cart_p)); /* @todo don't calculate */"
		"		vec3 sphe_p= vec3(dist, acos(cart_p.z/dist), atan(cart_p.y, cart_p.x));"
		"		vec2 s= sample(sphe_p);"
		"		intensity= max(0, intensity + (s.x - s.y*intensity)*dl/2.0);"
		"	}"
		"	intensity += rand(v_pos.xy + vec2(u_time/100.0, 0))*0.01;"
		"	gl_FragColor= vec4(color*intensity, 1.0);"
		"}"
		"\n";

	const std::size_t buf_size= 2048;
	char buf[buf_size];
	std::snprintf(buf,
		buf_size,
		"#version 120\n"
		"#define SAMPLE_COUNT %i\n"
		"#define WAVEFUNC(r, theta, phi, cos_theta, sin_theta) (%s)\n",
		sample_count,
		hydrogen_str.str);
	const GLchar* fs_src[]= { buf, fs_template_src };
	const GLsizei fs_src_count= sizeof(fs_src)/sizeof(*fs_src);

	Program::VolumeShader shd;
	createGlShaderProgram(shd.prog, shd.vs, shd.fs, 1, &vs_src, fs_src_count, fs_src);
	shd.timeLoc= glGetUniformLocation(shd.prog, "u_time");
	shd.phaseLoc= glGetUniformLocation(shd.prog, "u_phase");
	shd.colorLoc= glGetUniformLocation(shd.prog, "u_color");
	shd.transformLoc= glGetUniformLocation(shd.prog, "u_transform");
	return shd;
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
		font.whiteTexelUv= cast<Vec2f>(g_font.whitePixel)/font_size;

		// Calculate uv-coordinates for characters
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

		// Create vertically flipped RBGA texture from luminance data
		// This yields an OpenGL texture with conventional origin
		const Vec2i size= g_font.size;
		unsigned char rgba_data[size.x*size.y*4];
		for (int y= 0; y < size.y; ++y) {
			for (int x= 0; x < size.x; ++x) {
				int rgba_i= y*size.x*4 + x*4;
				int data_i= (size.y - y - 1)*size.x + x;
				rgba_data[rgba_i + 0]= 255;
				rgba_data[rgba_i + 1]= 255;
				rgba_data[rgba_i + 2]= 255;
				rgba_data[rgba_i + 3]= g_font.data[data_i];
			}
		}
		glGenTextures(1, &font.texId);
		glBindTexture(GL_TEXTURE_2D, font.texId);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(	GL_TEXTURE_2D, 0, GL_RGBA,
						size.x, size.y,
						0, GL_RGBA, GL_UNSIGNED_BYTE,
						rgba_data);
	}

	{ // Gui shader
		const GLchar* vs_src=
			"#version 120\n"
			"attribute vec2 a_pos;"
			"attribute vec2 a_uv;"
			"varying vec2 v_uv;"
			"void main() {"
			"	v_uv= a_uv;"
			"	gl_Position= vec4(a_pos, 0.0, 1.0);"
			"}\n";
		const GLchar* fs_src=
			"#version 120\n"
			"uniform sampler2D u_tex;"
			"uniform vec4 u_color;"
			"varying vec2 v_uv;"
			"void main() { gl_FragColor= texture2D(u_tex, v_uv)*u_color; }\n";

		Program::GuiShader& shd= prog.guiShader;
		createGlShaderProgram(shd.prog, shd.vs, shd.fs, 1, &vs_src, 1, &fs_src);
		shd.texLoc= glGetUniformLocation(shd.prog, "u_tex");
		shd.colorLoc= glGetUniformLocation(shd.prog, "u_color");
	}

	{ // Vbo used at rendering quads
		Program::QuadVbo& vbo= prog.vbo;
		glGenBuffers(1, &vbo.vboId);
		glBindBuffer(GL_ARRAY_BUFFER, vbo.vboId);
		glBufferData(GL_ARRAY_BUFFER, sizeof(Vec2f)*(4 + 4), NULL, GL_DYNAMIC_DRAW);
		glEnableVertexAttribArray(0); // Position
		glVertexAttribPointer(	0, 2, GL_FLOAT, GL_FALSE,
								sizeof(Vec2f)*2, BUFFER_OFFSET(0));
		glEnableVertexAttribArray(1); // Uv
		glVertexAttribPointer(	1, 2, GL_FLOAT, GL_FALSE,
								sizeof(Vec2f)*2, BUFFER_OFFSET(sizeof(Vec2f)));
	}

	{ // Program state
		prog.time= 0.0;
		prog.phase= 0.0;
		prog.sampleCount= 40;
		prog.resoMul= 0.5;
		prog.filtering= 0.0;
		prog.r= 0.4;
		prog.g= 0.8;
		prog.b= 1.0;
		prog.scale= 0.5;
		prog.n= 1;
		prog.l= 0;
		prog.m= 0;

		prog.shader= createVolumeShader(prog.sampleCount, prog.scale, prog.n, prog.l, prog.m);
		prog.fbo= createFbo(env.winSize*prog.resoMul, prog.filtering > 0.5);
	}

	{ // Setup initial GL state
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBindBuffer(GL_ARRAY_BUFFER, prog.vbo.vboId);
	}
}

void quit(Env& env, Program& prog)
{
	destroyFbo(prog.fbo);
	destroyGlShaderProgram(	prog.shader.prog,
							prog.shader.vs,
							prog.shader.fs);
	destroyGlShaderProgram(	prog.guiShader.prog,
							prog.guiShader.vs,
							prog.guiShader.fs);

	{ // Vbo
		glDeleteBuffers(1, &prog.vbo.vboId);
	}

	{ // Font
		glDeleteTextures(1, &prog.font.texId);
	}

	envQuit(env);
}

/// @note Uses currently bound vbo
void drawRect(	Vec2f ll,					Vec2f tr,
				Vec2f uv_ll= Vec2f(0, 0),	Vec2f uv_tr= Vec2f(1, 1))
{
	Vec2f v[4 + 4]= {
		ll, uv_ll,
		Vec2f(tr.x, ll.y), Vec2f(uv_tr.x, uv_ll.y),
		Vec2f(ll.x, tr.y), Vec2f(uv_ll.x, uv_tr.y),
		tr, uv_tr
	};
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(v), v);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

void frame(const Env& env, Program& prog)
{
	prog.time += env.dt;
	prog.phase += env.dt;

	struct Slider {
		const char* title;
		const float& min;
		const float& max;
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
			float v= clamp<float>((1.0 + x)/width()*(max - min) + min, min, max);
			return round(v, decimals);
		}
	};

	static Slider sliders[]= {
		{ "Phase",		0.0,	5.0,	prog.phase,			3, false },
		{ "Samples",	5,		150,	prog.sampleCount,	0, true },
		{ "Resolution",	0.01,	1.0,	prog.resoMul,		2, false },
		{ "Filtering",	0,		1,		prog.filtering,		0, false },
		{ "R",			0.0,	2.0,	prog.r,				3, false },
		{ "G",			0.0,	2.0,	prog.g,				3, false },
		{ "B",			0.0,	2.0,	prog.b,				3, false },
		{ "scale",		0.001,	1.0,	prog.scale,			3, true },
		{ "n",			1,		20,		prog.n,				0, true },
		{ "l",			0,		20,		prog.l,				0, true },
		{ "m",			0,		20,		prog.m,				0, true },
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
					destroyGlShaderProgram(	prog.shader.prog,
											prog.shader.vs,
											prog.shader.fs);
					prog.shader=
						createVolumeShader(prog.sampleCount, prog.scale, prog.n, prog.l, prog.m);
				}
			} else {
				slider_hover[i]= false;
			}
		}

		if (env.lmbDown && !slider_activity) {
			Vec2f smooth_delta= prev_delta*0.5 + (env.cursorDelta)*0.5;
			prev_delta= smooth_delta;
			
			rot += smooth_delta*2;
			rot.y= clamp<float>(rot.y, -tau/4, tau/4);
		}
	}

	// Slider texts
	const std::size_t slider_text_buf_size= 128;
	char slider_text[slider_count][slider_text_buf_size]= {};
	for (std::size_t i= 0; i < slider_count; ++i) {
		std::snprintf(	slider_text[i], slider_text_buf_size,
						"%s - %.*f", sliders[i].title, sliders[i].decimals, sliders[i].value);
	}

	{ // Adjust FBO
		Vec2i volume_reso= cast<Vec2i>(cast<Vec2f>(env.winSize)*prog.resoMul);
		bool volume_filtering= prog.filtering > 0.5;
		if (volume_reso != prog.fbo.reso || volume_filtering != prog.fbo.filtering) {
			destroyFbo(prog.fbo);
			prog.fbo= createFbo(volume_reso, volume_filtering);
		}
	}

	glClear(GL_COLOR_BUFFER_BIT);

	{ // Draw volume	
		// Draw to fbo
		glBindFramebuffer(GL_FRAMEBUFFER, prog.fbo.fboId);
		glViewport(0, 0, prog.fbo.reso.x, prog.fbo.reso.y);

		float s1= sin(rot.x), s2= sin(rot.y);
		float c1= cos(rot.x), c2= cos(rot.y);
		// Turntable-style rotation
		float transform[16]= {
			c1,			0,		s1,		0,
			-s1*s2,		c2,		c1*s2,	0,
			-c2*s1,		-s2,	c2*c1,	0,
			-c2*s1,		-s2,	c2*c1,	1 // Translation around origin
		};

		Program::VolumeShader& shd= prog.shader;
		glUseProgram(shd.prog);
		glUniform1f(shd.timeLoc, prog.time);
		glUniform1f(shd.phaseLoc, prog.phase);
		glUniform3f(shd.colorLoc, prog.r, prog.g, prog.b);
		glUniformMatrix4fv(shd.transformLoc, 1, GL_FALSE, transform);
		drawRect(Vec2f(-1, -1), Vec2f(1, 1));

		// Draw scaled fbo texture
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, env.winSize.x, env.winSize.y);
		glBindTexture(GL_TEXTURE_2D, prog.fbo.texId);
		glUseProgram(prog.guiShader.prog);
		glUniform1i(prog.guiShader.texLoc, 0);
		glUniform4f(prog.guiShader.colorLoc, 1.0, 1.0, 1.0, 1.0);
		drawRect(Vec2f(-1, -1), Vec2f(1, 1));
	}

	{ // Draw gui
		const Vec2f white_uv= prog.font.whiteTexelUv;
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glViewport(0, 0, env.winSize.x, env.winSize.y);
		glUseProgram(prog.guiShader.prog);
		glBindTexture(GL_TEXTURE_2D, prog.font.texId);
		glUniform1i(prog.guiShader.texLoc, 0);
		glUniform4f(prog.guiShader.colorLoc, 0.1, 0.1, 0.1, 0.3);

		// Background
		drawRect(	Vec2f(-1.0, 1.0 - slider_count*Slider::height()),
					Vec2f(Slider::width() - 1.0, 1.0),
					white_uv, white_uv);

		// Slider backgrounds
		for (std::size_t i= 0; i < slider_count; ++i) {
			Slider& s= sliders[i];
			float top= s.top(i);
			float bottom= s.bottom(i);
			float width= s.width()*s.fraction();

			if (!slider_hover[i])
				glUniform4f(prog.guiShader.colorLoc, 0.3, 0.3, 0.3, 0.6);
			else
				glUniform4f(prog.guiShader.colorLoc, 0.5, 0.5, 0.5, 0.8);

			drawRect(	Vec2f(-1.0, bottom), Vec2f(-1.0 + width, top),
						white_uv, white_uv);
		}

		// Slider texts
		glUniform4f(prog.guiShader.colorLoc, 0.8, 0.8, 0.8, 1.0);
		for (std::size_t s_i= 0; s_i < slider_count; ++s_i) {
			Slider& s= sliders[s_i];
			Vec2f pos= fitToGrid(Vec2f(-0.98, s.bottom(s_i)), env.winSize);

			for (std::size_t c_i= 0; c_i < std::strlen(slider_text[s_i]); ++c_i) {
				unsigned char ch= slider_text[s_i][c_i];
				Vec2f ch_size= cast<Vec2f>(g_font.charSize)/cast<Vec2f>(env.winSize)*2.0;
				Vec2f ch_pos(pos.x + ch_size.x*c_i, pos.y);

				Vec2f ch_tr= ch_pos + ch_size;
				Vec2f ll= prog.font.uv[ch];
				Vec2f tr= ll + prog.font.charUvSize;
				drawRect(ch_pos, ch_tr, ll, tr);
			}
		}
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

