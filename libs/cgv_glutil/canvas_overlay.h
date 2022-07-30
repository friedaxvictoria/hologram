#pragma once

#include <cgv_glutil/frame_buffer_container.h>
#include <cgv_glutil/overlay.h>
#include <cgv_glutil/2d/canvas.h>
#include <cgv_glutil/2d/shape2d_styles.h>
#include <cgv/gui/theme_info.h>

#include "lib_begin.h"

namespace cgv {
namespace glutil {

class CGV_API canvas_overlay : public overlay {
protected:
	int last_theme_idx = -1;

	frame_buffer_container fbc;

	canvas content_canvas, overlay_canvas;

	bool update_layout = false;
	bool has_damage = true;
	bool blend_overlay = false;
	GLboolean blending_was_enabled = false;

	void init_overlay_style(cgv::render::context& ctx);

	bool ensure_layout(cgv::render::context& ctx);

	bool ensure_theme();

	void begin_content(cgv::render::context& ctx, bool clear_frame_buffer = true);
	
	void end_content(cgv::render::context& ctx, bool keep_damage = false);

	void enable_blending();

	void disable_blending();

public:
	/// creates an overlay in the bottom left corner with zero size using a canvas for 2d drawing
	canvas_overlay();

	void clear(cgv::render::context& ctx);

	bool init(cgv::render::context& ctx);
	void draw(cgv::render::context& ctx);
	virtual void draw_content(cgv::render::context& ctx) = 0;

	void register_shader(const std::string& name, const std::string& filename);
};

typedef cgv::data::ref_ptr<canvas_overlay> canvas_overlay_ptr;

}
}

#include <cgv/config/lib_end.h>
