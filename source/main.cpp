// See unity.cpp for build instructions

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "env.hpp"
#include "fontdata.hpp"
#include "gl.hpp"

#define local_persist static

namespace qm {

const float sliderHeight= 0.05;
const float sliderWidth= 0.65;
struct Slider {
	const char* title;
	float min;
	float max;
	float* value;
	int decimals;
	bool recompile;

	static float top(std::size_t i) { return 1.0 - i*sliderHeight; }
	static float bottom(std::size_t i) { return top(i + 1); }
	bool pointInside(std::size_t i, Vec2f p) const
	{
		return	p.x  >= -1.0	&& p.x < sliderWidth - 1.0 &&
				p.y > bottom(i)	&& p.y < top(i);
	}
	float fraction() const
	{
		float f= (*value - min)/(max - min);
		return f < 1 ? f : 1;
	}
	float coordToValue(float x) const
	{
		float v= CLAMP((1.0 + x)/sliderWidth*(max - min) + min, min, max);
		return round(v, decimals);
	}
};

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

const std::size_t Program_maxSliders= 32;
const std::size_t Program_maxWaves= 2;
struct Program {
	VolumeShader shader;
	VolumeFbo fbo;
	Font font;
	GuiShader guiShader;
	QuadVbo vbo;
	float time;

	// Slider settings
	struct Wave {
		float amplitude;
		float phase; // Additional factor: e^(i*phase)
		float n, l, m; // Quantum numbers
	};
	float phase;
	float sampleCount;
	float resoMul; 
	float filtering; // bool
	float r, g, b;
	float complexColor; // bool
	float absorption;
	float cutoff;
	float scale;
	StackArray<Wave, Program_maxWaves> waves;
	StackArray<Slider, Program_maxSliders> sliders;
};

VolumeShader createVolumeShader(
		const int sample_count,
		const bool complex_color,
		const float absorption,
		const float cutoff,
		const float scale,
		const Program::Wave* waves,
		const std::size_t wave_count)
{
	/// @todo Remove
	testMath();

	assert(wave_count > 0);

	double bohr_radius= std::exp(scale*0.25) - 0.99999;
	const int poly_term_count= 30;
	// Formula for hydrogen wave function with parameters r, theta, and phi
	StackString<1024> hydrogen_amplitudes[Program_maxWaves]= {}; // Real multiplier
	StackString<64> hydrogen_phases[Program_maxWaves]= {}; // Complex phase
	for (std::size_t wave_i= 0; wave_i < wave_count; ++wave_i) {
		const Program::Wave& wave= waves[wave_i];
		StackString<1024>& hydrogen_amplitude= hydrogen_amplitudes[wave_i];
		StackString<64>& hydrogen_phase= hydrogen_phases[wave_i];
		float amplitude_adjust= wave.amplitude;
		int n= wave.n;
		int l= wave.l;
		int m= wave.m;
		assert(n > 0 && l >= 0);
		if (l > n - 1)
			l= n - 1;
		if (m > 0 && m > l)
			m= l;
		if (m < 0 && m < -l)
			m= -l;

		// Form a hydrogen wave function |nlm> in four parts
		// psi_nlm(r, theta, phi) = C*E*L*Y, where
		//   C = normalization factor sqrt[(2/(n*a_0))^3*(n - l - 1)!/(2n(n + l)!)]
		//   E = e^(-rho/l)*rho^l
		//   L = Generalized Laguerre Polynomial L(n - l - 1, 2l + 1, rho)
		//   Y = Spherical harmonic function Y(l, m, theta, phi) (phase is separated to different string)
		//   rho = 2r/(n*a_0)
		const std::size_t rho_str_size= 16;
		char rho_str[rho_str_size];
		std::snprintf(rho_str, rho_str_size, "%e*r", 2.0/bohr_radius/n);

		// C
		append(hydrogen_amplitude, "%e*%e",
			std::sqrt(	std::pow(2.0/(n*bohr_radius), 3)*
						fact(n - l - 1)/(2*n*fact(n + l))),
			amplitude_adjust);
		append(hydrogen_amplitude, "*");

		// E
		append(hydrogen_amplitude, "exp(-%s/2.0)*pow(%s, %i.0)", rho_str, rho_str, l);
		append(hydrogen_amplitude, "*");

		{ // L
			append(hydrogen_amplitude, "(");
			double laguerre_coeff[poly_term_count]= {};
			assert(n - l < poly_term_count);

			laguerre(laguerre_coeff, n - l - 1, 2*l + 1);
			for (int i= 0; i < poly_term_count; ++i) {
				if (laguerre_coeff[i] == 0.0)
					continue;

				append(hydrogen_amplitude, "+%e", laguerre_coeff[i]);
				if (i != 0)
					append(hydrogen_amplitude, "*pow(%s, %i.0)", rho_str, i);
			}
			append(hydrogen_amplitude, ")");
		}
		append(hydrogen_amplitude, "*");

		{ // Y
			if (m != 0)
				append(hydrogen_amplitude,
					"%s*pow(abs(sin_theta), %i.0)*",
					(std::abs(m) % 2 ? "sign(sin_theta)" : "1.0"),
					std::abs(m));
			append(hydrogen_amplitude, "(");

			double spherical_coeff[poly_term_count]= {};
			assert(l - 1 < poly_term_count);
			sphericalHarmonics(spherical_coeff, l, m);
			for (int i= 0; i < poly_term_count; ++i) {
				if (spherical_coeff[i] == 0.0)
					continue;

				append(hydrogen_amplitude, "+(%e)", spherical_coeff[i]);
				if (i != 0)
					append(hydrogen_amplitude,
						"*%s*pow(abs(cos_theta), %i.0)",
						(i % 2 ? "sign(cos_theta)" : "1.0"),
						i);
			}
			append(hydrogen_amplitude, ")");

			append(hydrogen_phase, "%i.0*phi + (%e)", m, wave.phase);
		}

		//std::printf("Wave function:\n%s\n", hydrogen_amplitude.str);
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
		"void main()"
		"{"
		"	vec3 n= normalize(v_normal);"
		"	vec3 color= u_color;"
		"	vec3 intensity= vec3(0.0, 0.0, 0.0);"
		"	float last_P= 0.0;"
		"	const float dl= 2.0/SAMPLE_COUNT;"
		"	for (int i= 0; i < SAMPLE_COUNT; ++i) {"
		"		vec3 cart_p= v_pos + n*2.0*float(SAMPLE_COUNT - i - 1)/float(SAMPLE_COUNT);"
		"		float r= sqrt(dot(cart_p, cart_p));"
		"		float phi= atan(cart_p.y, cart_p.x);"
		"		float cos_theta= cart_p.z/r;"
		"		float theta= acos(cos_theta);"
		"		float sin_theta= sin(theta);"
		"		float total_real= 0;"
		"		float total_imag= 0;"
		"		CALC_TOTAL_WAVEFUNC;"
		"		float total_amplitude= total_real*total_real + total_imag*total_imag;"
		"		float total_complex_phase= atan(total_imag, total_real);"
		"		float P= total_amplitude*total_amplitude;"
		"		if (P < CUTOFF) P= 0;"
		"\n#if COMPLEX_COLOR == 1\n"
		"		vec3 emission= P*normalize(vec3(0.5 - 0.5*cos(total_complex_phase), 0.2, 0.5 + 0.5*sin(total_complex_phase)));"
		"\n#else\n"
		"		vec3 emission= P*color;"
		"\n#endif\n"
		"		float absorption= P*ABSORPTION_MUL;"
		"		intensity=	intensity + (emission - absorption*intensity)*dl/2.0;"
		"		intensity= max(vec3(0.0, 0.0, 0.0), intensity);"
		"		last_P= P;"
		"	}"
		"	intensity += vec3(1.0, 1.0, 1.0)*rand(v_pos.xy*u_time)*0.015;"
		"	gl_FragColor= vec4(intensity, 1.0);"
		"}"
		"\n";

	StackString<1024*Program_maxWaves> calc_total_wavefunc_define= {};
	append(calc_total_wavefunc_define, "#define CALC_TOTAL_WAVEFUNC ");
	for (int i= 0; i < (int)wave_count; ++i) {
		if (waves[i].amplitude <= 0.001)
			continue;
		append(calc_total_wavefunc_define,
			"float a_%i= 2.3456*(%s);" // Can be negative
			"float p_%i= (%s);" // Not taking account possible negative amplitude
			"total_real += a_%i*cos(p_%i);"
			"total_imag += a_%i*sin(p_%i);",
			i, hydrogen_amplitudes[i].str, i, hydrogen_phases[i].str, i, i, i, i, i);
	}

	StackString<1024*Program_maxWaves + 256> buf= {};
	append(buf,
		"#version 120\n"
		"#define SAMPLE_COUNT %i\n"
		"#define COMPLEX_COLOR %i\n"
		"#define ABSORPTION_MUL %e\n"
		"#define CUTOFF %e\n"
		"%s\n",
		sample_count,
		complex_color,
		absorption,
		cutoff,
		calc_total_wavefunc_define.str);
	const GLchar* fs_src[]= { buf.str, fs_template_src };
	const GLsizei fs_src_count= sizeof(fs_src)/sizeof(*fs_src);

	VolumeShader shd;
	createGlShaderProgram(shd.prog, shd.vs, shd.fs, 1, &vs_src, fs_src_count, fs_src);
	shd.timeLoc= glGetUniformLocation(shd.prog, "u_time");
	shd.phaseLoc= glGetUniformLocation(shd.prog, "u_phase");
	shd.colorLoc= glGetUniformLocation(shd.prog, "u_color");
	shd.transformLoc= glGetUniformLocation(shd.prog, "u_transform");
	return shd;
}

VolumeFbo createFbo(Vec2i reso, bool filtering)
{
	GLenum filter= filtering ? GL_LINEAR : GL_NEAREST;

	VolumeFbo fbo;
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

void destroyFbo(VolumeFbo& fbo)
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glDeleteFramebuffers(1, &fbo.fboId);
	glDeleteTextures(1, &fbo.texId);
}

void addWave(Program& prog)
{
	Program::Wave w= {};
	w.n= 1;
	if (prog.waves.size == 0)
		w.amplitude= 0.5;
	push(prog.waves, w);

	Slider amplitude= { "Amplitude", 0.0, 1.0, &last(prog.waves).amplitude, 3, true };
	Slider phase= { "Complex phase", 0.0, (float)tau, &last(prog.waves).phase, 3, true };
	Slider n= { "n", 1, 15, &last(prog.waves).n, 0, true };
	Slider l= { "l", 0, 14, &last(prog.waves).l, 0, true };
	Slider m= { "m", -14, 14, &last(prog.waves).m, 0, true };
	push(prog.sliders, amplitude);
	push(prog.sliders, phase);
	push(prog.sliders, n);
	push(prog.sliders, l);
	push(prog.sliders, m);
}

void removeWave(Program& prog)
{
	pop(prog.waves);
	pop(prog.sliders); // n
	pop(prog.sliders); // l
	pop(prog.sliders); // m
}

void init(Env& env, Program& prog)
{
	env= envInit();
	queryGlFuncs();

	{ // Font
		Font& font= prog.font;
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

		GuiShader& shd= prog.guiShader;
		createGlShaderProgram(shd.prog, shd.vs, shd.fs, 1, &vs_src, 1, &fs_src);
		shd.texLoc= glGetUniformLocation(shd.prog, "u_tex");
		shd.colorLoc= glGetUniformLocation(shd.prog, "u_color");
	}

	{ // Vbo used at rendering quads
		QuadVbo& vbo= prog.vbo;
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
		prog.r= 1.0;
		prog.g= 0.6;
		prog.b= 0.4;
		prog.complexColor= 0.0;
		prog.absorption= 0.0;
		prog.cutoff= 0.0;
		prog.scale= 0.5;

		Slider default_sliders[] = {
			{ "Time",			0.0,	5.0,	&prog.phase,			3, false },
			{ "Samples",		5,		150,	&prog.sampleCount,		0, true },
			{ "Resolution",		0.01,	1.0,	&prog.resoMul,			2, false },
			{ "Filtering",		0,		1,		&prog.filtering,		0, false },
			{ "R",				0.0,	2.0,	&prog.r,				3, false },
			{ "G",				0.0,	2.0,	&prog.g,				3, false },
			{ "B",				0.0,	2.0,	&prog.b,				3, false },
			{ "Complex color",	0,		1,		&prog.complexColor,		0, true },
			{ "Absorption",		0.0,	1.0,	&prog.absorption,		3, true },
			{ "Cutoff",			0.0,	3.0,	&prog.cutoff,			3, true },
			{ "Scale",			0.001,	1.0,	&prog.scale,			3, true },
		};
		const std::size_t default_slider_count= sizeof(default_sliders)/sizeof(*default_sliders);
		for (std::size_t i= 0; i < default_slider_count; ++i)
			push(prog.sliders, default_sliders[i]);

		addWave(prog);
		addWave(prog);

		prog.shader=
			createVolumeShader(
				prog.sampleCount,
				prog.complexColor,
				prog.absorption,
				prog.cutoff,
				prog.scale, 
				prog.waves.data, prog.waves.size);
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

	bool slider_hover[Program_maxSliders]= {};

	local_persist Vec2f prev_delta;
	local_persist Vec2f rot;

	{ // User interaction
		bool slider_activity= false;
		for (std::size_t i= 0; i < prog.sliders.size; ++i) {
			Slider& s= prog.sliders.data[i];
			if (s.pointInside(i, env.anchorPos)) {
				bool value_changed= false;
				if (env.lmbDown) {
					float new_value= s.coordToValue(env.cursorPos.x);
					value_changed= new_value != *s.value;
					*s.value= new_value;
				}
				slider_hover[i]= true;
				slider_activity= true;

				if (value_changed && s.recompile) {
					destroyGlShaderProgram(	prog.shader.prog,
											prog.shader.vs,
											prog.shader.fs);
					prog.shader=
						createVolumeShader(
							prog.sampleCount,
							prog.complexColor,
							prog.absorption,
							prog.cutoff,
							prog.scale, 
							prog.waves.data, prog.waves.size);
				}
			} else {
				slider_hover[i]= false;
			}
		}

		if (env.lmbDown && !slider_activity) {
			Vec2f smooth_delta= prev_delta*0.5 + (env.cursorDelta)*0.5;
			prev_delta= smooth_delta;
			
			rot += smooth_delta*2;
			rot.y= CLAMP(rot.y, -tau/4, tau/4);
		}
	}

	// Slider texts
	const std::size_t slider_text_buf_size= 128;
	char slider_text[Program_maxSliders][slider_text_buf_size]= {};
	for (std::size_t i= 0; i < prog.sliders.size; ++i) {
		Slider& slider= prog.sliders.data[i];
		std::snprintf(	slider_text[i], slider_text_buf_size,
						"%s - %.*f",
						slider.title,
						slider.decimals,
						*slider.value);
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

		VolumeShader& shd= prog.shader;
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
		drawRect(	Vec2f(-1.0, 1.0 - prog.sliders.size*sliderHeight),
					Vec2f(sliderWidth - 1.0, 1.0),
					white_uv, white_uv);

		// Slider backgrounds
		for (std::size_t i= 0; i < prog.sliders.size; ++i) {
			Slider& s= prog.sliders.data[i];
			float top= s.top(i);
			float bottom= s.bottom(i);
			float width= sliderWidth*s.fraction();

			if (!slider_hover[i])
				glUniform4f(prog.guiShader.colorLoc, 0.3, 0.3, 0.3, 0.6);
			else
				glUniform4f(prog.guiShader.colorLoc, 0.5, 0.5, 0.5, 0.8);

			drawRect(	Vec2f(-1.0, bottom), Vec2f(-1.0 + width, top),
						white_uv, white_uv);
		}

		// Slider texts
		glUniform4f(prog.guiShader.colorLoc, 0.8, 0.8, 0.8, 1.0);
		for (std::size_t s_i= 0; s_i < prog.sliders.size; ++s_i) {
			Slider& s= prog.sliders.data[s_i];
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
	qm::Env env= {};
	qm::Program prog= {};
	qm::init(env, prog);

	while (!env.quitRequested) {
		qm::envUpdate(env);
		qm::frame(env, prog);
	}

	qm::quit(env, prog);
}

