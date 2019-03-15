#include "plot_base.h"
#include <cgv/render/shader_program.h>
#include <cgv/signal/rebind.h>

namespace cgv {
	namespace plot {

std::vector<const char*> plot_base::font_names;
std::string plot_base::font_name_enum_def;

///
plot_base::tick_batch_info::tick_batch_info(int _ai, int _aj, bool _primary, unsigned _first_vertex, unsigned _first_label) 
	: ai(_ai), aj(_aj), primary(_primary), first_vertex(_first_vertex), first_label(_first_label)
{
}

/// set tick config defaults
tick_config::tick_config(bool primary)
{
	step = primary ? 5.0f : 1.0f;
	type = TT_LINE;
	line_width = primary ? 1.5f : 1.0f;
	length = primary ? 2.0f : 1.0f;
	label = primary;
	precision = -1;
}

axis_config::axis_config() : primary_ticks(true), secondary_ticks(false), color(0.3f,0.3f,0.3f)
{
	log_scale = false;
	line_width = 3.0f;
}

domain_config::domain_config(unsigned nr_axes) : color(0.85f,0.85f,0.85f), axis_configs(nr_axes)
{
	show_domain = true;
	fill = true;
	label_font_index = -1;
	label_font_size = 16.0f;
	label_ffa = cgv::media::font::FFA_BOLD_ITALIC;
}

plot_base_config::plot_base_config(const std::string& _name) : name(_name)
{
	show_plot = true;

	point_size = 3;
	point_color = rgb(1,0,0);

	stick_width = 2;
	stick_color = rgb(1,1,0);

	bar_percentual_width = 0.75f;
	bar_outline_width = 1;
	bar_color = rgb(0,1,1);
	bar_outline_color = rgb(1,1,1);

	configure_chart(CT_BAR_CHART);
}

/// configure the sub plot to a specific chart type
void plot_base_config::configure_chart(ChartType chart_type)
{
	switch (chart_type) {
	case CT_POINT : 
	case CT_LINE_CHART:
		show_points = true;
		show_sticks = false;
		show_bars = false;
		break;
	case CT_BAR_CHART:
		show_points = false;
		show_sticks = false;
		show_bars = true;
		break;
	}
}


plot_base_config::~plot_base_config()
{
}

/// check whether tick information has to be updated
bool plot_base::tick_render_information_outofdate() const
{
	if (last_dom_cfg.axis_configs.size() != get_domain_config_ptr()->axis_configs.size())
		return true;
	for (unsigned ai = 0; ai < last_dom_cfg.axis_configs.size(); ++ai) {
		const axis_config& ac = last_dom_cfg.axis_configs[ai];
		const axis_config& ao = get_domain_config_ptr()->axis_configs[ai];
		if (ac.log_scale != ao.log_scale)
			return true;
		for (unsigned ti = 0; ti < 2; ++ti) {
			const tick_config& tc = ti == 0 ? ac.primary_ticks : ac.secondary_ticks;
			const tick_config& to = ti == 0 ? ao.primary_ticks : ao.secondary_ticks;
			if (tc.length != to.length)
				return true;
			if (tc.step != to.step)
				return true;
			if (tc.type != to.type)
				return true;
			if (tc.label != to.label)
				return true;
			if (tc.precision != to.precision)
				return true;
		}
	}
	if (last_dom_cfg.label_font_size!= get_domain_config_ptr()->label_font_size)
		return true;
	if (last_dom_cfg.label_font_index != get_domain_config_ptr()->label_font_index)
		return true;
	if (last_dom_cfg.label_ffa != get_domain_config_ptr()->label_ffa)
		return true;
	return false;
}

/// ensure that tick render information is current
void plot_base::ensure_tick_render_information()
{
	if (tick_render_information_outofdate()) {
		tick_vertices.clear();
		tick_labels.clear();
		tick_batches.clear();
		compute_tick_render_information();
		last_dom_cfg = *get_domain_config_ptr();
	}
}

float log_conform_add(float v0, float v1, bool log_scale, float v_min, float v_max)
{
	if (!log_scale)
		return v0 + v1;
	float q = (v0 - 0.5f*(v_min+v_max)) / (v_max - v_min);
	return pow(10.0f, (q*(log(v_max) - log(v_min)) + 0.5f*(log(v_min) + log(v_max))) / log(10.0f));
}

void plot_base::collect_tick_geometry(int ai, int aj, const float* dom_min_pnt, const float* dom_max_pnt, const float* extent)
{
	axis_config& ac = get_domain_config_ptr()->axis_configs[ai];
	for (unsigned ti=0; ti<2; ++ti) {
		tick_config& tc = ti == 0 ? ac.primary_ticks : ac.secondary_ticks;
		if (tc.type == TT_NONE)
			return;
		tick_batches.push_back(tick_batch_info(ai, aj, ti == 0, tick_vertices.size(), tick_labels.size()));

		axis_config& ao = get_domain_config_ptr()->axis_configs[aj];
		// compute domain extent in both coordinate directions
		float dei = dom_max_pnt[ai] - dom_min_pnt[ai];
		float dej = dom_max_pnt[aj] - dom_min_pnt[aj];
		float dci = dom_min_pnt[ai] + 0.5f*dei;
		float dcj = dom_min_pnt[aj] + 0.5f*dej;

		float min_val = dom_min_pnt[ai];
		float max_val = dom_max_pnt[ai];
		if (ac.log_scale) {
			min_val = log10(min_val);
			max_val = log10(max_val);
		}
		int min_i = (int)((min_val - fmod(min_val, tc.step)) / tc.step);
		int max_i = (int)((max_val - fmod(max_val, tc.step)) / tc.step);

		// ignore ticks on domain boundary
		if (min_i * tc.step - min_val < std::numeric_limits<float>::epsilon())
			++min_i;
		if (max_i * tc.step - max_val > -std::numeric_limits<float>::epsilon())
			--max_i;

		float dash_length = tc.length*0.01f*dej;
		if (extent[ai] < extent[aj])
			dash_length *= extent[ai] / extent[aj];

		float s_min = dom_min_pnt[aj];
		float s_max = dom_max_pnt[aj];

		for (int i = min_i; i <= max_i; ++i) {
			vec2 c;
			c[ai] = (float)(i*tc.step);
			if (ac.log_scale)
				c[ai] = pow(10.0f, c[ai]);
			else // ignore ticks on axes
				if (fabs(c[ai]) < std::numeric_limits<float>::epsilon())
					continue;

			std::string label_str;
			if (tc.label)
				label_str = cgv::utils::to_string(c(ai));
			switch (tc.type) {
			case TT_DASH:
				// generate label
				if (!label_str.empty()) {
					// left label
					c[aj] = log_conform_add(s_min, -0.5f*dash_length, ao.log_scale, dom_min_pnt[aj], dom_max_pnt[aj]);
					tick_labels.push_back(label_info(c, label_str, ai == 0 ? cgv::render::TA_TOP : cgv::render::TA_RIGHT));
					// right label
					c[aj] = log_conform_add(s_max,  0.5f*dash_length, ao.log_scale, dom_min_pnt[aj], dom_max_pnt[aj]);
					tick_labels.push_back(label_info(c, label_str, ai == 0 ? cgv::render::TA_BOTTOM : cgv::render::TA_LEFT));
				}
				// left tick
				c[aj] = s_min;
				tick_vertices.push_back(c);
				c[aj] = log_conform_add(s_min,  dash_length, ao.log_scale, dom_min_pnt[aj], dom_max_pnt[aj]);
				tick_vertices.push_back(c);

				// right tick
				c[aj] = s_max;
				tick_vertices.push_back(c);
				c[aj] = log_conform_add(s_max, -dash_length, ao.log_scale, dom_min_pnt[aj], dom_max_pnt[aj]);
				tick_vertices.push_back(c);

				// non log axis tick
				if (!ao.log_scale && s_min + 0.5f*dash_length < 0 && s_max - 0.5f*dash_length > 0) {
					c[aj] = -0.5f*dash_length;
					tick_vertices.push_back(c);
					c[aj] = 0.5f*dash_length;
					tick_vertices.push_back(c);
				}
				break;
			case TT_LINE:
			case TT_PLANE:
				// generate label
				if (!label_str.empty()) {
					// left label
					c[aj] = log_conform_add(s_min, -0.5f*dash_length, ao.log_scale, dom_min_pnt[aj], dom_max_pnt[aj]);
					tick_labels.push_back(label_info(c, label_str, ai == 0 ? cgv::render::TA_TOP : cgv::render::TA_RIGHT));
					// right label
					c[aj] = log_conform_add(s_max, 0.5f*dash_length, ao.log_scale, dom_min_pnt[aj], dom_max_pnt[aj]);
					tick_labels.push_back(label_info(c, label_str, ai == 0 ? cgv::render::TA_BOTTOM : cgv::render::TA_LEFT));
				}
				c[aj] = s_min; tick_vertices.push_back(c);
				if (tc.label) {
					c(aj) -= 0.5f*dash_length;
					if (ao.log_scale) {
						float q = (c[aj] - dcj) / dej;
						c[aj] = pow(10.0f, (q*(log(dom_max_pnt[aj]) - log(dom_min_pnt[aj])) + 0.5f*(log(dom_min_pnt[aj]) + log(dom_max_pnt[aj]))) / log(10.0f));
					}
					tick_labels.push_back(label_info(c, label_str, ai == 0 ? cgv::render::TA_TOP : cgv::render::TA_RIGHT));
				}
				c[aj] = s_max; tick_vertices.push_back(c);
				break;
			}
		}
		tick_batches.back().vertex_count = tick_vertices.size() - tick_batches.back().first_vertex;
		tick_batches.back().label_count = tick_labels.size() - tick_batches.back().first_label;
	}
}


void plot_base::ensure_font_names()
{
	if (font_names.empty()) {
		cgv::media::font::enumerate_font_names(font_names);
		if (font_names.empty())
			return;
		font_name_enum_def = std::string("enums='") + font_names[0];
		for (unsigned i = 1; i < font_names.size(); ++i) {
			font_name_enum_def += ',';
			font_name_enum_def += font_names[i];
		}
		font_name_enum_def += "'";
	}
	if (!label_font) {
		get_domain_config_ptr()->label_font_index = 0;
		for (auto iter = font_names.begin(); iter != font_names.end(); ++iter)
			if (std::string(*iter) == "Times New Roman") {
				get_domain_config_ptr()->label_font_index = iter - font_names.begin();
				break;
			}
		on_font_selection();
	}
}

void plot_base::on_font_selection()
{
	label_font = cgv::media::font::find_font(font_names[get_domain_config_ptr()->label_font_index]);
	on_font_face_selection();
}

void plot_base::on_font_face_selection()
{
	label_font_face = label_font->get_font_face(get_domain_config_ptr()->label_ffa);
}

void plot_base::set_uniforms(cgv::render::context& ctx, cgv::render::shader_program& prog, unsigned i)
{
	prog.set_uniform(ctx, "point_color", ref_sub_plot_config(i).point_color);
	prog.set_uniform(ctx, "stick_color", ref_sub_plot_config(i).stick_color);
	prog.set_uniform(ctx, "bar_color", ref_sub_plot_config(i).bar_color);
	prog.set_uniform(ctx, "bar_outline_color", ref_sub_plot_config(i).bar_outline_color);
}

plot_base::plot_base(unsigned nr_axes) : dom_cfg(nr_axes), last_dom_cfg(0)
{
	dom_cfg_ptr = &dom_cfg;
}

/// configure the label font
void plot_base::set_label_font(float font_size, cgv::media::font::FontFaceAttributes ffa, const std::string& font_name)
{
	get_domain_config_ptr()->label_font_size = font_size;
	get_domain_config_ptr()->label_ffa = ffa;
	if (!font_name.empty()) {
		for (auto iter = font_names.begin(); iter != font_names.end(); ++iter)
			if (font_name == *iter) {
				get_domain_config_ptr()->label_font_index = iter - font_names.begin();
				on_font_selection();
				break;
			}
	}
	else
		on_font_face_selection();
}

/// return const pointer to domain configuration
const domain_config* plot_base::get_domain_config_ptr() const
{
	return dom_cfg_ptr;
}

/// return pointer to domain configuration
domain_config* plot_base::get_domain_config_ptr()
{
	return dom_cfg_ptr;
}

/// set the domain configuration to an external configuration in order to synch several plots, if set to null, the internal domain config is used again
void plot_base::set_domain_config_ptr(domain_config* _new_ptr)
{
	if (_new_ptr)
		dom_cfg_ptr = _new_ptr;
	else
		dom_cfg_ptr = &dom_cfg;
}


unsigned plot_base::get_nr_sub_plots() const
{
	return configs.size();
}

plot_base_config& plot_base::ref_sub_plot_config(unsigned i)
{
	return *configs[i];
}

/// ensure tick computation
void plot_base::init_frame(cgv::render::context& ctx)
{
	ensure_tick_render_information();
}

void plot_base::create_plot_gui(cgv::base::base* bp, cgv::gui::provider& p)
{
	const char* axis_names = "xyz";

	ensure_font_names();

	bool open = p.begin_tree_node("domain", get_domain_config_ptr()->show_domain, false, "level=3;w=70;align=' '");
	p.add_member_control(bp, "show", get_domain_config_ptr()->show_domain, "toggle", "w=40", " ");
	p.add_member_control(bp, "fill", get_domain_config_ptr()->fill, "toggle", "w=40");
	if (open) {
		p.align("\a");
		p.add_member_control(bp, "fill color", get_domain_config_ptr()->color);
		for (unsigned i = 0; i < get_domain_config_ptr()->axis_configs.size(); ++i) {
			axis_config& ac = get_domain_config_ptr()->axis_configs[i];
			bool show = p.begin_tree_node(std::string(1, axis_names[i]) + " axis", ac.color, false, "level=3;w=100;align=' '");
			p.add_member_control(bp, "log", ac.log_scale, "toggle", "w=50");
			if (show) {
				p.align("\a");
				p.add_member_control(bp, "width", ac.line_width, "value_slider", "min=1;max=20;log=true;ticks=true");
				p.add_member_control(bp, "color", ac.color);
				char* tn[2] = { "primary tick", "secondary tick" };
				tick_config* tc[2] = { &ac.primary_ticks, &ac.secondary_ticks };
				for (unsigned ti = 0; ti < 2; ++ti) {
					bool vis = p.begin_tree_node(tn[ti], tc[ti]->label, false, "level=3;w=100;align=' '");
					p.add_member_control(bp, "label", tc[ti]->label, "toggle", "w=60");
					if (vis) {
						p.align("\a");
						p.add_member_control(bp, "type", tc[ti]->type, "dropdown", "enums='none,dash,line,plane'");
						p.add_member_control(bp, "step", tc[ti]->step, "value");
						p.add_member_control(bp, "width", tc[ti]->line_width, "value_slider", "min=1;max=20;log=true;ticks=true");
						p.add_member_control(bp, "length", tc[ti]->length, "value_slider", "min=1;max=20;log=true;ticks=true");
						p.add_member_control(bp, "precision", tc[ti]->precision, "value_slider", "min=-1;max=5;ticks=true");
						p.align("\b");
						p.end_tree_node(tc[ti]->label);
					}
				}
				p.align("\b");
				p.end_tree_node(ac.color);
			}
		}
		p.add_member_control(bp, "font_size", get_domain_config_ptr()->label_font_size, "value_slider", "min=8;max=40;log=false;ticks=true");
		p.add_member_control(bp, "font", (cgv::type::DummyEnum&)get_domain_config_ptr()->label_font_index, "dropdown", font_name_enum_def);
		connect_copy(p.find_control((cgv::type::DummyEnum&)get_domain_config_ptr()->label_font_index)->value_change,
			cgv::signal::rebind(this, &plot_base::on_font_selection));
		p.add_member_control(bp, "face", get_domain_config_ptr()->label_ffa, "dropdown", "enums='normal,bold,italics,bold italics'");
		connect_copy(p.find_control(get_domain_config_ptr()->label_ffa)->value_change,
			cgv::signal::rebind(this, &plot_base::on_font_face_selection));
		p.align("\b");
		p.end_tree_node(get_domain_config_ptr()->show_domain);
	}
}

void plot_base::create_config_gui(cgv::base::base* bp, cgv::gui::provider& p, unsigned i)
{
	plot_base_config& pbc = ref_sub_plot_config(i);
	p.add_member_control(bp, "name", pbc.name);
	bool show = p.begin_tree_node("points", pbc.show_points, false, "level=3;w=100;align=' '");
	p.add_member_control(bp, "show", pbc.show_points, "toggle", "w=50");
	if (show) {
		p.align("\a");
			p.add_member_control(bp, "size", pbc.point_size, "value_slider", "min=1;max=20;log=true;ticks=true");
			p.add_member_control(bp, "color", pbc.point_color);
		p.align("\b");
		p.end_tree_node(pbc.show_points);
	}
	show = p.begin_tree_node("sticks", pbc.show_sticks, false, "level=3;w=100;align=' '");
	p.add_member_control(bp, "show", pbc.show_sticks, "toggle", "w=50");
	if (show) {
		p.align("\a");
			p.add_member_control(bp, "width", pbc.stick_width, "value_slider", "min=1;max=20;log=true;ticks=true");
			p.add_member_control(bp, "color", pbc.stick_color);
		p.align("\b");
		p.end_tree_node(pbc.show_sticks);
	}
	show = p.begin_tree_node("bars", pbc.show_bars, false, "level=3;w=100;align=' '");
	p.add_member_control(bp, "show", pbc.show_bars, "toggle", "w=50");
	if (show) {
		p.align("\a");
			p.add_member_control(bp, "width", pbc.bar_percentual_width, "value_slider", "min=0.01;max=1;log=true;ticks=true");
			p.add_member_control(bp, "fill", pbc.bar_color);
			p.add_member_control(bp, "outline_width", pbc.bar_outline_width, "value_slider", "min=0;max=20;log=true;ticks=true");
			p.add_member_control(bp, "outline", pbc.bar_outline_color);
		p.align("\b");
		p.end_tree_node(pbc.show_bars);
	}
}

void plot_base::create_gui(cgv::base::base* bp, cgv::gui::provider& p)
{
	create_plot_gui(bp, p);
	for (unsigned i=0; i<get_nr_sub_plots(); ++i) {
		plot_base_config& pbc = ref_sub_plot_config(i);
		bool show = p.begin_tree_node(std::string("configure ")+pbc.name, pbc.name, false, "level=3;w=100;align=' '");
		p.add_member_control(bp, "show", pbc.show_plot, "toggle", "w=50");
		if (show) {
			p.align("\a");
				create_config_gui(bp, p, i);
			p.align("\b");
			p.end_tree_node(pbc.name);
		}
	}
}


	}
}

#ifdef REGISTER_SHADER_FILES
#include <cgv/base/register.h>
#include <plot_shader_inc.h>
#endif
