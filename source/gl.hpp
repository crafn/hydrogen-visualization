#ifndef QM_GL_HPP
#define QM_GL_HPP

#include <GL/gl.h>
#include "env.hpp"

// Partial OpenGL 2.1 interface

typedef char GLchar;

#define GL_FRAGMENT_SHADER 0x8B30
#define GL_VERTEX_SHADER 0x8B31
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82

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
typedef void (*GlUniform3f)(GLuint, GLfloat, GLfloat, GLfloat);
GlUniform3f glUniform3f;

namespace qm {

inline
void queryGlFuncs()
{
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
	glUniform3f= (GlUniform3f)queryGlFunc("glUniform3f");
}

inline
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

inline
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

} // qm

#endif // QM_GL_HPP
