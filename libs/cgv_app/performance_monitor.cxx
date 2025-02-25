#include "performance_monitor.h"

#include <cgv/gui/theme_info.h>
#include <cgv/math/ftransform.h>

namespace cgv {
namespace app {

performance_monitor::performance_monitor() {

	set_name("Performance Monitor");
	block_events = false;
	blend_overlay = true;
	gui_options.allow_stretch = false;

	layout.padding = 13; // 10px plus 3px border
	layout.total_size = ivec2(180, 90);

	set_overlay_alignment(AO_START, AO_END);
	set_overlay_stretch(SO_NONE);
	set_overlay_margin(ivec2(-3));
	set_overlay_size(layout.total_size);

	register_shader("rectangle", cgv::g2d::canvas::shaders_2d::rectangle);
	register_shader("line", cgv::g2d::canvas::shaders_2d::line);

	bar_renderer = cgv::g2d::generic_2d_renderer(cgv::g2d::canvas::shaders_2d::rectangle);
}

void performance_monitor::clear(cgv::render::context& ctx) {

	canvas_overlay::clear(ctx);

	cgv::g2d::ref_msdf_font(ctx, -1);
	cgv::g2d::ref_msdf_gl_canvas_font_renderer(ctx, -1);

	label_font.destruct(ctx);

	bar_renderer.destruct(ctx);
	bars.destruct(ctx);
}

bool performance_monitor::self_reflect(cgv::reflect::reflection_handler& _rh) {

	return false;
}

bool performance_monitor::handle_event(cgv::gui::event& e) {

	return false;
}

void performance_monitor::on_set(void* member_ptr) {

	if(member_ptr == &layout.total_size[0] || member_ptr == &layout.total_size[1]) {
		
	}

	if(member_ptr == &show_plot) {
		layout.total_size.y() = show_plot ? 90 : 55;
		set_overlay_size(layout.total_size);
	}

	if(member_ptr == &show_background || member_ptr == &invert_color) {
		auto ctx_ptr = get_context();
		if(ctx_ptr)
			init_styles(*ctx_ptr);
	}

	if(member_ptr == &monitor.enabled) {
		if(monitor.enabled) {
			monitor.timer.restart();
			monitor.total_frame_count = 0u;
			monitor.interval_frame_count = 0u;
			monitor.last_seconds_since_start = 0.0;
			monitor.running_time = 0.0;
		}
	}

	update_member(member_ptr);
	post_damage();
}

bool performance_monitor::init(cgv::render::context& ctx) {

	bool success = canvas_overlay::init(ctx);

	success &= bar_renderer.init(ctx);

	if(success)
		init_styles(ctx);

	cgv::g2d::msdf_font& font = cgv::g2d::ref_msdf_font(ctx, 1);
	cgv::g2d::ref_msdf_gl_canvas_font_renderer(ctx, 1);

	if(font.is_initialized()) {
		texts.set_msdf_font(&font);
		texts.set_font_size(text_font_size);
	}

	label_font.set_font_face(cgv::g2d::msdf_font::FF_LIGHT);
	label_font.init(ctx);

	if(label_font.is_initialized()) {
		labels.set_msdf_font(&label_font);
		labels.set_font_size(label_font_size);
	}

	plot_color_map.add_color_point(0.0f, rgb(0.5f, 1.0f, 0.5f));
	plot_color_map.add_color_point(0.25f, rgb(0.0f, 0.9f, 0.0f));
	plot_color_map.add_color_point(0.5f, rgb(0.8f, 0.9f, 0.0f));
	plot_color_map.add_color_point(1.0f, rgb(0.9f, 0.0f, 0.0f));

	return success;
}

void performance_monitor::init_frame(cgv::render::context& ctx) {

	if(ensure_layout(ctx)) {
		ivec2 container_size = get_overlay_size();
		layout.update(container_size);
		create_texts();
		create_labels();
	}

	if(ensure_theme())
		init_styles(ctx);

	if(monitor.enabled) {
		if(show_plot)
			update_plot();
		update_stats_texts();
	}
}

void performance_monitor::draw_content(cgv::render::context& ctx) {

	begin_content(ctx);

	ivec2 container_size = get_overlay_size();

	auto& font_renderer = cgv::g2d::ref_msdf_gl_canvas_font_renderer(ctx);
	auto& rect_prog = content_canvas.enable_shader(ctx, "rectangle");

	// draw container background
	if(show_background) {
		container_style.apply(ctx, rect_prog);
		content_canvas.draw_shape(ctx, ivec2(0), container_size);
	}

	if(show_plot) {
		// draw plot border
		border_style.apply(ctx, rect_prog);
		content_canvas.draw_shape(ctx, layout.plot_rect.pos() - 1, layout.plot_rect.size() + 2);

		// draw plot bars
		bar_renderer.render(ctx, content_canvas, cgv::render::PT_POINTS, bars, bar_style);

		// draw line
		auto& line_prog = content_canvas.enable_shader(ctx, "line");

		const auto& r = layout.plot_rect;
		ivec2 a(r.x() + 12, r.box.get_center().y());
		ivec2 b = a;
		b.x() = r.x1();

		line_style.apply(ctx, line_prog);
		content_canvas.draw_shape2(ctx, a, b);
		content_canvas.disable_current_shader(ctx);

		// draw labels
		font_renderer.render(ctx, content_canvas, labels, label_style);
	}

	// draw text
	font_renderer.render(ctx, content_canvas, texts, text_style);

	end_content(ctx);
}

void performance_monitor::after_finish(cgv::render::context& ctx) {

	if(monitor.enabled) {
		++monitor.total_frame_count;
		++monitor.interval_frame_count;
		
		double seconds_since_start = monitor.timer.get_elapsed_time();
		monitor.delta_time = seconds_since_start - monitor.last_seconds_since_start;
		
		monitor.running_time += monitor.delta_time;

		monitor.last_seconds_since_start = seconds_since_start;

		if(monitor.running_time >= monitor.interval) {
			monitor.avg_fps = (double)monitor.interval_frame_count / monitor.running_time;
			monitor.running_time = 0.0;
			monitor.interval_frame_count = 0u;
		}
	}
}

void performance_monitor::enable_monitoring(bool enabled) {
	monitor.enabled = enabled;
	on_set(&monitor.enabled);
}

void performance_monitor::create_gui_impl() {

	add_member_control(this, "Enable", monitor.enabled, "check", "w=110", " ");
	add_member_control(this, "Show Plot", show_plot, "check", "w=78");
	add_member_control(this, "Measure Interval (s)", monitor.interval, "value_slider", "min=0.01;max=1;step=0.01;ticks=true");
	
	add_member_control(this, "Background", show_background, "check", "w=100", " ");
	add_member_control(this, "Invert Color", invert_color, "check", "w=88");
}

void performance_monitor::init_styles(cgv::render::context& ctx) {
	// get theme colors
	auto& ti = cgv::gui::theme_info::instance();
	rgba background_color = rgba(ti.background(), 1.0f);
	rgba group_color = rgba(ti.group(), 1.0f);
	rgb border_color = ti.text();

	if(invert_color) {
		border_color.R() = pow(1.0f - pow(border_color.R(), 2.2f), 1.0f/2.2f);
		border_color.G() = pow(1.0f - pow(border_color.G(), 2.2f), 1.0f/2.2f);
		border_color.B() = pow(1.0f - pow(border_color.B(), 2.2f), 1.0f/2.2f);
	}

	// configure style for the container rectangle
	container_style.fill_color = group_color;
	container_style.border_color = background_color;
	container_style.border_width = 3.0f;
	container_style.feather_width = 0.0f;

	// configure style for the border rectangle
	border_style = container_style;
	border_style.fill_color = show_background ? rgba(ti.text_background(), 1.0f) : rgba(0.0f);
	border_style.border_color = rgba(border_color, 1.0);
	border_style.border_width = 1.0f;
	border_style.feather_width = 0.0f;
	border_style.use_blending = true;

	line_style.use_blending = true;
	line_style.fill_color = rgba(border_color, invert_color ? 0.666f : 0.333f);
	line_style.feather_width = 0.0f;
	line_style.dash_length = 10.0f;
	
	bar_style.use_fill_color = false;
	bar_style.feather_width = 0.0f;

	// configure text style
	float label_border_alpha = 0.0f;
	float border_width = 0.25f;
	if(!show_background) {
		label_border_alpha = 1.0f;
		border_width = 0.0f;
	}

	text_style.fill_color = rgba(border_color, 1.0f);
	text_style.border_color = rgba(border_color, label_border_alpha);
	text_style.border_width = border_width;
	text_style.feather_origin = 0.5f;
	text_style.use_blending = true;

	label_style = text_style;
	label_style.feather_width = 0.5f;
}

void performance_monitor::create_texts() {

	texts.clear();

	const int line_spacing = static_cast<int>(1.25f* text_font_size);

	ivec2 caret_pos = ivec2(layout.content_rect.x(), layout.content_rect.y1() - (int)text_font_size);
	texts.add_text("Frames per Second:", caret_pos, cgv::render::TA_BOTTOM_LEFT);
	caret_pos.y() -= line_spacing;
	texts.add_text("Frametime (ms):", caret_pos, cgv::render::TA_BOTTOM_LEFT);
	
	caret_pos = ivec2(layout.content_rect.x1(), layout.content_rect.y1() - (int)text_font_size);
	texts.add_text("", caret_pos, cgv::render::TA_BOTTOM_RIGHT);
	caret_pos.y() -= line_spacing;
	texts.add_text("", caret_pos, cgv::render::TA_BOTTOM_RIGHT);
}

void performance_monitor::update_stats_texts() {

	if(texts.size() > 3) {
		std::stringstream ss;
		ss.precision(2);
		ss << std::fixed;
		ss << monitor.avg_fps;

		std::string str = ss.str();
		texts.set_text(2, ss.str());

		ss.str(std::string());
		if(monitor.avg_fps < 0.001f)
			ss << "-";
		else
			ss << 1000.0 / monitor.avg_fps;

		texts.set_text(3, ss.str());

		post_damage();
	}
}

void performance_monitor::create_labels() {
	
	labels.clear();

	ivec2 caret_pos = ivec2(layout.plot_rect.x(), layout.plot_rect.y1() + 2);
	labels.add_text("30", caret_pos, cgv::render::TA_TOP_LEFT);
	caret_pos.y() = layout.plot_rect.box.get_center().y() + 1;
	labels.add_text("60", caret_pos, cgv::render::TA_LEFT);
	caret_pos.y() = layout.plot_rect.y();
	labels.add_text("120", caret_pos, cgv::render::TA_BOTTOM_LEFT);
}

void performance_monitor::update_plot() {

	ivec2 plot_size = layout.plot_rect.size();

	float a = static_cast<float>(1000.0 * monitor.delta_time / 33.333333333);
	float b = std::min(a, 1.0f);
	float bar_height = plot_size.y() * b;
	bar_height = std::max(bar_height, 1.0f);

	rgb bar_color = a > 1.0f ? rgb(0.7f, 0.0f, 0.0f) : plot_color_map.interpolate_color(b);

	if(bars.get_render_count() < plot_size.x()) {
		for(auto& position : bars.position)
			position.x() -= 1.0f;

		int x = layout.plot_rect.x1() - 1;
		x = std::max(x, 0);
		float bar_x = static_cast<float>(x);
		bars.add(vec2(bar_x, static_cast<float>(layout.plot_rect.y())), vec2(1.0f, bar_height), bar_color);

	} else {
		for(size_t i = 0; i < bars.position.size() - 1; ++i) {
			bars.size[i].y() = bars.size[i + 1].y();
			bars.color[i] = bars.color[i + 1];
		}

		bars.size.back() = vec2(1.0f, bar_height);
		bars.color.back() = bar_color;
	}

	bars.set_out_of_date();
	post_damage();
}

}
}
