#pragma once

#include <cgv/gui/event_handler.h>
#include <cgv/gui/provider.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv/render/drawable.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/texture.h>
#include <cgv/render/vertex_buffer.h>
#include <cgv_gl/box_renderer.h>
#include <cgv_glutil/frame_buffer_container.h>
#include <cgv_glutil/overlay.h>
#include <cgv_glutil/shader_library.h>

#include "transfer_function.h"
#include "2d/draggables_collection.h"

#include "lib_begin.h"

namespace cgv {
namespace glutil{

class CGV_API transfer_function_editor : public overlay {
protected:
	bool show;
	
	std::string file_name;
	std::string save_file_name;
	bool has_unsaved_changes = false;

	bool mouse_is_on_overlay;
	bool show_cursor;
	ivec2 cursor_pos;
	std::string cursor_drawtext;
	cgv::media::font::font_face_ptr cursor_font_face;

	cgv::glutil::frame_buffer_container fbc;
	cgv::glutil::shader_library shaders;

	float opacity_scale_exponent;
	cgv::type::DummyEnum resolution;

	bool show_histogram;
	rgba histogram_color;
	rgba histogram_border_color;
	unsigned histogram_border_width;
	float histogram_smoothing;

	struct layout_attributes {
		int padding;
		int total_height;
		int color_scale_height;

		// dependent members
		rect editor_rect;
		rect color_scale_rect;

		void update(const ivec2& parent_size) {

			editor_rect.set_pos(ivec2(padding) + ivec2(0, color_scale_height + 1)); // add one for a small border
			editor_rect.set_size(parent_size - 2 * padding - ivec2(0, color_scale_height + 1));

			color_scale_rect.set_pos(ivec2(padding));
			color_scale_rect.set_size(ivec2(parent_size.x() - 2*padding, color_scale_height));
		}
	} layout;
	
	struct point : public draggable {
		vec2 val;
		rgb col;

		point() {
			size = vec2(8.0f);
			position_is_center = true;
			constraint_reference = CR_CENTER;
		}

		void update_val(const layout_attributes& la, const float scale_exponent) {

			vec2 p = pos - la.editor_rect.pos();
			val = p / la.editor_rect.size();

			val = cgv::math::clamp(val, 0.0f, 1.0f);
			val.y() = cgv::math::clamp(std::pow(val.y(), scale_exponent), 0.0f, 1.0f);
		}

		void update_pos(const layout_attributes& la, const float scale_exponent) {

			val = cgv::math::clamp(val, 0.0f, 1.0f);

			vec2 t = val;

			t.y() = cgv::math::clamp(std::pow(t.y(), 1.0f / scale_exponent), 0.0f, 1.0f);

			pos = la.editor_rect.pos() + t * la.editor_rect.size();
		}

		bool is_inside(const vec2& mp) const {

			float dist = length(mp - center());
			return dist <= size.x();
		}

		ivec2 get_render_position() const {
			return ivec2(pos + 0.5f);
		}

		ivec2 get_render_size() const {
			return 2 * ivec2(size);
		}
	};

	texture bg_tex;
	texture tf_tex;

	template<typename PosType, typename ColType>
	struct plain_geometry {
		struct vertex_type {
			PosType pos;
			ColType col;
		};

		std::vector<vertex_type> vertices;

		type_descriptor pos_type_descriptor = cgv::render::element_descriptor_traits<PosType>::get_type_descriptor(PosType());
		type_descriptor col_type_descriptor = cgv::render::element_descriptor_traits<ColType>::get_type_descriptor(ColType());

		vertex_buffer vb;
		attribute_array_binding aab;

		size_t size() { return vertices.size(); }

		void clear(context& ctx) {
			vertices.clear();

			if(vb.is_created())
				vb.destruct(ctx);
			if(aab.is_created())
				aab.destruct(ctx);
		}

		bool create(context& ctx, const shader_program& prog) {
			bool success = true;
			success &= vb.create(ctx, &(vertices[0]), vertices.size());
			success &= aab.create(ctx);
			success &= aab.set_attribute_array(ctx, prog.get_position_index(), pos_type_descriptor, vb, 0, vertices.size(), sizeof(vertex_type));
			success &= aab.set_attribute_array(ctx, prog.get_color_index(), col_type_descriptor, vb, sizeof(PosType), vertices.size(), sizeof(vertex_type));
			return success;
		}

		void add(const PosType& pos, const ColType& col) {
			vertices.push_back({ pos, col });
		}

		void render(context& ctx, PrimitiveType type, shader_program& prog) {
			render(ctx, type, 0, size(), prog);
		}

		void render(context& ctx, PrimitiveType type, int offset, size_t count, shader_program& prog) {
			if(aab.is_created()) {
				prog.enable(ctx);
				aab.enable(ctx);
				GLenum mode = gl::map_to_gl(type);
				glDrawArrays(mode, (GLint)offset, (GLsizei)count);
				aab.disable(ctx);
				prog.disable(ctx);
			}
		}
	};

	struct tf_container {
		cgv::glutil::draggables_collection<point> points;

		transfer_function tf;
		texture tex;

		texture hist_tex;
		unsigned hist_max;

		plain_geometry<vec2, rgba> triangles;
		plain_geometry<vec2, rgb> lines;

		void reset() {
			points.clear();

			point p;
			p.val = vec2(0.0f);
			p.col = rgb(0.0f);
			points.add(p);

			p = point();
			p.val = vec2(1.0f);
			p.col = rgb(1.0f);
			points.add(p);
		}
	} tfc;

	void add_point(const vec2& pos);
	void remove_point(const point* ptr);
	point* get_hit_point(const vec2& pos);
	
	void init_transfer_function_texture(context& ctx);

	void handle_drag();
	void handle_drag_end();
	void sort_points();
	void update_point_positions();
	void update_transfer_function(bool is_data_change);
	bool update_geometry();

	bool load_from_xml(const std::string& file_name);
	bool save_to_xml(const std::string& file_name);

public:
	transfer_function_editor();
	std::string get_type_name() const { return "transfer_function_editor"; }

	bool on_exit_request();
	void clear(cgv::render::context& ctx);

	bool self_reflect(cgv::reflect::reflection_handler& _rh);
	void stream_help(std::ostream& os) {}

	bool handle_event(cgv::gui::event& e);
	void on_set(void* member_ptr);

	bool init(cgv::render::context& ctx);
	void init_frame(cgv::render::context& ctx);
	void draw(cgv::render::context& ctx);
	
	void create_gui();
	void create_gui(cgv::gui::provider& p);

	void is_visible(bool visible);
	void toggle_visibility();

	texture& ref_tex();

	bool set_histogram(const std::vector<unsigned>& data);
};

}
}

#include <cgv/config/lib_end.h>
