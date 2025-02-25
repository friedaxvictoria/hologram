#pragma once

#include "msdf_gl_font_renderer.h"
#include "canvas.h"

#include "lib_begin.h"

namespace cgv {
namespace g2d {

class CGV_API msdf_gl_canvas_font_renderer;

extern CGV_API msdf_gl_canvas_font_renderer& ref_msdf_gl_canvas_font_renderer(cgv::render::context& ctx, int ref_count_change = 0);

class CGV_API msdf_gl_canvas_font_renderer : public msdf_gl_font_renderer {
public:
	bool enable(cgv::render::context& ctx, canvas& cvs, msdf_text_geometry& tg, const shape2d_style& style) {
		bool res = msdf_gl_font_renderer::enable(ctx, cvs.get_resolution(), tg, style);
		if(res)
			cvs.set_view(ctx, prog);
		return res;
	}

	void draw(cgv::render::context& ctx, msdf_text_geometry& tg, size_t offset = 0, int count = -1);

	void draw(cgv::render::context& ctx, canvas& cvs, msdf_text_geometry& tg, size_t offset = 0, int count = -1);

	bool render(cgv::render::context& ctx, canvas& cvs, msdf_text_geometry& tg, const shape2d_style& style, size_t offset = 0, int count = -1);
};

}
}

#include <cgv/config/lib_end.h>
