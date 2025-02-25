#include "msdf_text_geometry.h"

namespace cgv {
namespace g2d {

msdf_text_geometry::msdf_text_geometry() {
	msdf_font_ptr = nullptr;
	state_out_of_date = true;

	render_font_size = 32.0f;
};

msdf_text_geometry::~msdf_text_geometry() {
	clear();
};

void msdf_text_geometry::clear() {
	state_out_of_date = true;

	texts.clear();
	vertices.clear();

	if(ssbo != 0) {
		glDeleteBuffers(1, &ssbo);
		ssbo = 0;
	}
}

void msdf_text_geometry::set_msdf_font(msdf_font* ptr, bool update_texts) {
	msdf_font_ptr = ptr;

	if(msdf_font_ptr && update_texts) {
		for(size_t i = 0; i < texts.size(); ++i)
			set_text(i, texts[i].str);
	}
}

void msdf_text_geometry::set_text(unsigned i, const std::string& text) {
	if(i < texts.size()) {
		texts[i].str = text;
		texts[i].size.x() = compute_length(text);

		update_offsets_and_counts();

		state_out_of_date = true;
	}
}

void msdf_text_geometry::set_position(unsigned i, const ivec2& position) {
	if(i < texts.size())
		texts[i].position = position;
}

void msdf_text_geometry::set_alignment(unsigned i, const cgv::render::TextAlignment alignment) {
	if(i < texts.size())
		texts[i].alignment = alignment;
}

void msdf_text_geometry::set_angle(unsigned i, const float angle) {
	if(i < texts.size())
		texts[i].angle = angle;
}

msdf_text_geometry::vec2 msdf_text_geometry::get_text_render_size(unsigned i) const {

	if(i < texts.size())
		return render_font_size * texts[i].size;

	return vec2(0.0f);
}

void msdf_text_geometry::add_text(const std::string& str, const ivec2& position, const cgv::render::TextAlignment alignment, float angle) {
	text_info text;
	text.str = str;
	text.position = position;
	text.size = vec2(compute_length(text.str), 1.0f);
	text.alignment = alignment;
	text.angle = angle;

	end_text(text);
	state_out_of_date = true;
}

bool msdf_text_geometry::create(cgv::render::context& ctx) {
	create_vertex_data();

	if(ssbo != 0) {
		glDeleteBuffers(1, &ssbo);
		ssbo = 0;
	}

	glGenBuffers(1, &ssbo);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, vertices.size() * sizeof(vertex_type), vertices.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	state_out_of_date = false;
	return true;
}

bool msdf_text_geometry::enable(cgv::render::context& ctx) {
	if(state_out_of_date) {
		create(ctx);
	}

	if(!msdf_font_ptr) {
		std::cerr << "msdf_text_geometry::enable: msdf_font is not specified" << std::endl;
		return false;
	}

	bool success = msdf_font_ptr->enable(ctx);
	if(success)
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);
	return success;
}

void msdf_text_geometry::disable(cgv::render::context& ctx) {
	msdf_font_ptr->disable(ctx);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

float msdf_text_geometry::compute_length(const std::string& str) const {
	float length = 0.0f;
	float acc_advance = 0.0f;

	for(char c : str) {
		const msdf_font::glyph_info& g = msdf_font_ptr->get_glyph_info(static_cast<unsigned char>(c));
		length = acc_advance + g.size.x();
		acc_advance += g.advance;
	}

	return length;
}

void msdf_text_geometry::end_text(text_info text) {
	text.offset = 0;
	if(texts.size() > 0) {
		const text_info& last_text = *texts.rbegin();
		text.offset = int(last_text.offset + last_text.str.size());
	}
	texts.push_back(text);
}

void msdf_text_geometry::update_offsets_and_counts() {
	int offset = 0;

	for(text_info& text : texts) {
		text.offset = offset;
		offset += int(text.str.size());
	}
}

void msdf_text_geometry::add_vertex(const vec4& pos, const vec4& txc) {
	vertices.push_back({ pos, txc });
}

void msdf_text_geometry::create_vertex_data() {
	vertices.clear();

	if(!msdf_font_ptr)
		return;

	for(text_info& text : texts) {
		float acc_advance = 0.0f;

		for(char c : text.str) {
			const msdf_font::glyph_info& g = msdf_font_ptr->get_glyph_info(static_cast<unsigned char>(c));

			vec2 position = g.position + vec2(acc_advance, 0.0f);
			vec2 size = g.size;
			acc_advance += g.advance;

			add_vertex(vec4(position.x(), position.y(), size.x(), size.y()), g.texcoords);
		}
	}
}

}
}
