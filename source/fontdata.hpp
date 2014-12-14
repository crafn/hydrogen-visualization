#ifndef QM_FONTDATA_HPP
#define QM_FONTDATA_HPP

#include "util.hpp"

namespace qm {

struct FontData {
	const char* chars;
	Vec2i charSize;
	Vec2i size;
	const unsigned char* data;
	Vec2i whitePixel;
};

extern const FontData g_font;

} // qm

#endif // QM_FONTDATA_HPP
