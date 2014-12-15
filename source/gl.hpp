#ifndef QM_GL_HPP
#define QM_GL_HPP

#include <GL/gl.h>
#include "env.hpp"

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

// Required GL 2.1 features

#define GL_ARRAY_BUFFER 0x8892
#define GL_DYNAMIC_DRAW 0x88E8
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_VERTEX_SHADER 0x8B31
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82

typedef char GLchar;
typedef intptr_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

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
typedef GLint (*GlGetUniformLocation)(GLuint, const GLchar*);
GlGetUniformLocation glGetUniformLocation;
typedef void (*GlUniform1f)(GLuint, GLfloat);
GlUniform1f glUniform1f;
typedef void (*GlUniform3f)(GLuint, GLfloat, GLfloat, GLfloat);
GlUniform3f glUniform3f;
typedef void (*GlUniform4f)(GLuint, GLfloat, GLfloat, GLfloat, GLfloat);
GlUniform4f glUniform4f;
typedef void (*GlUniformMatrix4fv)(GLint, GLsizei, GLboolean, const GLfloat*);
GlUniformMatrix4fv glUniformMatrix4fv;
typedef void (*GlUniform1i)(GLint, GLint);
GlUniform1i glUniform1i;
typedef void (*GlGenBuffers)(GLsizei, GLuint*);
GlGenBuffers glGenBuffers;
typedef void (*GlBindBuffer)(GLenum, GLuint);
GlBindBuffer glBindBuffer;
typedef void (*GlBufferData)(GLenum, GLsizeiptr, const GLvoid*, GLenum);
GlBufferData glBufferData;
typedef void (*GlBufferSubData)(GLenum, GLintptr, GLsizeiptr, const GLvoid*);
GlBufferSubData glBufferSubData;
typedef void (*GlDeleteBuffers)(GLsizei, const GLuint*);
GlDeleteBuffers glDeleteBuffers;
typedef void (*GlEnableVertexAttribArray)(GLuint);
GlEnableVertexAttribArray glEnableVertexAttribArray;
typedef void (*GlVertexAttribPointer)(GLuint, GLint, GLenum, GLboolean, GLsizei, const GLvoid*);
GlVertexAttribPointer glVertexAttribPointer;

// Required GL 3 features

#define GL_FRAMEBUFFER 0x8D40
#define GL_COLOR_ATTACHMENT0 0x8CE0

typedef void (*GlGenFramebuffers)(GLsizei, GLuint*);
GlGenFramebuffers glGenFramebuffers;
typedef void (*GlBindFramebuffer)(GLenum, GLuint);
GlBindFramebuffer glBindFramebuffer;
typedef void (*GlFramebufferTexture2D)(GLenum, GLenum, GLenum, GLuint, GLint);
GlFramebufferTexture2D glFramebufferTexture2D;
typedef void (*GlDeleteFramebuffers)(GLsizei, GLuint*);
GlDeleteFramebuffers glDeleteFramebuffers;

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
	glUniform4f= (GlUniform4f)queryGlFunc("glUniform4f");
	glUniformMatrix4fv= (GlUniformMatrix4fv)queryGlFunc("glUniformMatrix4fv");
	glUniform1i= (GlUniform1i)queryGlFunc("glUniform1i");
	glGenBuffers= (GlGenBuffers)queryGlFunc("glGenBuffers");
	glBindBuffer= (GlBindBuffer)queryGlFunc("glBindBuffer");
	glBufferData= (GlBufferData)queryGlFunc("glBufferData");
	glBufferSubData= (GlBufferSubData)queryGlFunc("glBufferSubData");
	glDeleteBuffers= (GlDeleteBuffers)queryGlFunc("glDeleteBuffers");
	glEnableVertexAttribArray= (GlEnableVertexAttribArray)queryGlFunc("glEnableVertexAttribArray");
	glVertexAttribPointer= (GlVertexAttribPointer)queryGlFunc("glVertexAttribPointer");

	glGenFramebuffers= (GlGenFramebuffers)queryGlFunc("glGenFramebuffers");
	glBindFramebuffer= (GlBindFramebuffer)queryGlFunc("glBindFramebuffer");
	glFramebufferTexture2D= (GlFramebufferTexture2D)queryGlFunc("glFramebufferTexture2D");
	glDeleteFramebuffers= (GlDeleteFramebuffers)queryGlFunc("glDeleteFramebuffers");
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

inline
void checkGlErrors(const char* tag)
{
	GLenum error= glGetError();
	if (error!= GL_NO_ERROR)
		std::printf("GL Error (%s): %i\n", tag, error);
}

inline
void createGlShaderProgram(	GLuint& prog, GLuint& vs, GLuint& fs,
							GLsizei vs_count, const GLchar** vs_src,
							GLsizei fs_count, const GLchar** fs_src)
{
	{ // Vertex
		vs= glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vs, vs_count, vs_src, NULL);
		glCompileShader(vs);
		checkShaderStatus(vs);
	}
	{ // Fragment
		fs= glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fs, fs_count, fs_src, NULL);
		glCompileShader(fs);
		checkShaderStatus(fs);
	}
	{ // Shader program
		prog= glCreateProgram();
		glAttachShader(prog, vs);
		glAttachShader(prog, fs);
		glLinkProgram(prog);
		checkProgramStatus(prog);
	}
}

inline
void destroyGlShaderProgram(GLuint prog, GLuint vs, GLuint fs)
{
	glDetachShader(prog, vs);
	glDeleteShader(vs);

	glDetachShader(prog, fs);
	glDeleteShader(fs);

	glDeleteProgram(prog);
}

} // qm

#endif // QM_GL_HPP
