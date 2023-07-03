#include "holo_view_interactor.h"
#include <cgv/gui/application.h>
#include <cgv/base/group.h>
#include <cgv/math/geom.h>
#include <cgv/math/ftransform.h>
#include <cgv/gui/dialog.h>
#include <libs/cg_gamepad/gamepad_server.h>
#include <cgv_reflect_types/math/fvec.h>
#include <cgv_reflect_types/media/color.h>
#include <cgv/render/shader_program.h>
#include <cgv/utils/scan.h>
#include <cgv/utils/scan_enum.h>
#include <cgv/utils/ostream_printf.h>
#include <cgv/gui/trigger.h>
#include <cgv/gui/animate.h>
#include <cgv/signal/rebind.h>
#include <cgv/gui/key_event.h>
#include <cgv/gui/mouse_event.h>
#include <cgv/media/image/image_writer.h>
#include <cgv/type/variant.h>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <filesystem>
#include <plugins/cg_fltk/fltk_gl_view.h>

using namespace cgv::math;
using namespace cgv::signal;
using namespace cgv::gui;
using namespace cgv::utils;
using namespace cgv::render;
using namespace cgv::render::gl;
using namespace cgv::base;

// enable (set to 1) in case of evaluation - hole visualisation is disabled and writing the evaluation's results to a .csv file is enabled
#define EVAL 0
// enable (set to 1) in case compute shader are evaluated - unnecessary code that is needed for vertex warping is not compiled
#define COMPUTE 0

#define SMP_ENUMS "bitmap,pixels,arrow"
#define HOLO_ENUMS "single,quilt,volume"

cgv::reflect::enum_reflection_traits<holo_view_interactor::StorageMode>
get_reflection_traits(const holo_view_interactor::StorageMode&)
{
	return cgv::reflect::enum_reflection_traits<holo_view_interactor::StorageMode>("single,quilt,volume");
}

void holo_view_interactor::set_default_values()
{
	stereo_view::set_default_values();
	eye_distance = 0.3f;
	two_d_enabled = false;
}

void holo_view_interactor::check_emulation_active()
{
	emulation_active = false;
	for (unsigned i = 0; i < 6; ++i) {
		if (plus_key_values[i] > 0 || plus_key_down[i] || minus_key_values[i] > 0 || minus_key_down[i]) {
			emulation_active = true;
			return;
		}
	}
}

void holo_view_interactor::plus_key_action(int i, cgv::gui::KeyAction action)
{
	switch (action) {
	case cgv::gui::KA_RELEASE:
		plus_key_toggle_time[i] = cgv::gui::trigger::get_current_time();
		plus_key_down[i] = false;
		break;
	case cgv::gui::KA_PRESS:
		plus_key_toggle_time[i] = cgv::gui::trigger::get_current_time();
		plus_key_down[i] = true;
		break;
	case cgv::gui::KA_REPEAT:
		plus_key_down[i] = true;
		break;
	}
	check_emulation_active();
}
void holo_view_interactor::minus_key_action(int i, cgv::gui::KeyAction action)
{
	switch (action) {
	case cgv::gui::KA_RELEASE:
		minus_key_toggle_time[i] = cgv::gui::trigger::get_current_time();
		minus_key_down[i] = false;
		break;
	case cgv::gui::KA_PRESS:
		minus_key_toggle_time[i] = cgv::gui::trigger::get_current_time();
		minus_key_down[i] = true;
		break;
	case cgv::gui::KA_REPEAT:
		minus_key_down[i] = true;
		break;
	}
	check_emulation_active();
}

void holo_view_interactor::timer_event(double t, double dt)
{
	if (emulation_active) {
		unsigned i;
		// update key values
		for (i = 0; i < 6; ++i) {
			if (plus_key_down[i] || plus_key_values[i] > 0) {
				double dt = t - plus_key_toggle_time[i];
				if (!plus_key_down[i])
					dt = -dt;
				plus_key_values[i] += (float)(0.2 * dt);
				if (plus_key_values[i] > 1)
					plus_key_values[i] = 1;
				if (plus_key_values[i] < 0)
					plus_key_values[i] = 0;
			}
			if (minus_key_down[i] || minus_key_values[i] > 0) {
				double dt = t - minus_key_toggle_time[i];
				if (!minus_key_down[i])
					dt = -dt;
				minus_key_values[i] += (float)(0.2 * dt);
				if (minus_key_values[i] > 1)
					minus_key_values[i] = 1;
				if (minus_key_values[i] < 0)
					minus_key_values[i] = 0;
			}
			emulation_axes[i] = plus_key_values[i] - minus_key_values[i];
		}
		left_stick[0] = emulation_axes[0];
		left_stick[1] = emulation_axes[1];
		right_stick[0] = emulation_axes[2];
		right_stick[1] = emulation_axes[3];
		check_emulation_active();
	}
	else {
		if (!(left_stick.length() > deadzone || right_stick.length() > deadzone ||
			  fabs(trigger[1] - trigger[0]) > deadzone))
			return;
	}
	dvec3 x, y, z;
	put_coordinate_system(x, y, z);
	if (left_stick.length() > deadzone) {
		float mode_sign = (right_mode == 1 ? -1.0f : 1.0f);
		switch (left_mode) {
		case 0:
			rotate(4 * mode_sign * dt / rotate_sensitivity * left_stick[1],
				   -4 * dt / rotate_sensitivity * left_stick[0], right_mode == 1 ? 0.0 : get_depth_of_focus());
			on_rotation_change();
			break;
		case 1:
			pan(-5 * mode_sign * dt * get_y_extent_at_focus() / pan_sensitivity * left_stick[0],
				-5 * mode_sign * dt * get_y_extent_at_focus() / pan_sensitivity * left_stick[1]);
			update_vec_member(view::focus);
			break;
		}
		post_redraw();
	}
	if (right_stick.length() > deadzone) {
		switch (right_mode) {
		case 0:
			if (fabs(right_stick[0]) > fabs(right_stick[1])) {
				roll(-4 * dt / rotate_sensitivity * right_stick[0]);
				on_rotation_change();
			}
			else {
				zoom(pow(2.0f, -10 * dt * right_stick[1] / zoom_sensitivity));
				update_member(&y_extent_at_focus);
			}
			break;
		case 1:
			if (fabs(right_stick[0]) > fabs(right_stick[1])) {
				set_y_view_angle(get_y_view_angle() + 200 * dt / zoom_sensitivity * right_stick[0]);
				if (get_y_view_angle() < 1)
					set_y_view_angle(1);
				if (get_y_view_angle() > 160)
					set_y_view_angle(160);
				update_member(&y_view_angle);
			}
			else {
				set_focus(get_focus() - 20 * dt * get_y_extent_at_focus() * right_stick[1] * z / zoom_sensitivity);
				update_vec_member(view::focus);
			}
			break;
		}
		post_redraw();
	}
}

///
holo_view_interactor::holo_view_interactor(const char* name)
	: node(name), quilt_depth_buffer("[D]"), volume_depth_buffer("[D]"), layered_depth_tex("flt32[D]"),
	  layered_color_tex("flt32[R,G,B,A]"), volume_holo_tex("flt32[R,G,B,A]"), quilt_holo_tex("flt32[R,G,B,A]")
{
	enable_messages = true;
	use_gamepad = true;
	gamepad_emulation = false;
	emulation_active = false;
	for (unsigned i = 0; i < 6; ++i) {
		emulation_axes[i] = 0;
		plus_key_values[i] = 0;
		plus_key_down[i] = false;
		plus_key_toggle_time[i] = 0;
		minus_key_values[i] = 0;
		minus_key_toggle_time[i] = 0;
		minus_key_down[i] = false;
	}
	last_do_viewport_splitting = do_viewport_splitting = false;
	viewport_shrinkage = 2;
	deadzone = 0.03f;
	gamepad_attached = false;
	left_mode = right_mode = 0;
	left_stick = right_stick = trigger = cgv::math::fvec<float, 2>(0.0f);
	connect(cgv::gui::get_animation_trigger().shoot, this, &holo_view_interactor::timer_event);

	fix_view_up_dir = false;
	stereo_translate_in_model_view = false;
	set_default_values();
	z_near_derived = z_near;
	z_far_derived = z_far;
	clip_relative_to_extent = false;
	show_focus = false;
	check_for_click = -1;
	pan_sensitivity = zoom_sensitivity = rotate_sensitivity = 1;
	last_x = last_y = -1;
}
/// return the type name
std::string holo_view_interactor::get_type_name() const { return "holo_view_interactor"; }

/// overload to stream help information to the given output stream
void holo_view_interactor::stream_help(std::ostream& os)
{
	os << "holo_view_interactor:\n\a"
	   << "stereo[F4], stereo mode[s-F4], z_near<c-N,a-N>, z_far<c-F,a-F>, select view dir{c-X|Y|Z,sc-X|Y|Z}\n"
	   << "view all{c-Spc}, show focus[s-F], set focus{click LMB}, pan<RMB>, zoom<MMB>, rotate<LMB>, roll<s-LMB>";
	if (fix_view_up_dir)
		os << " disabled";
	os << "\n";
	os << "zoom to focus<MW> or <PgUp,PgDn>, dolly zoom<s-MW>, eye separation<a+MW>, parallax zero plane<c-MW>";
	if (gamepad_attached) {
		os << "\nPAD: view all{start},focus{A}, select view dir{DPAD|Bumpers},navigate<sticks>";
	}
	os << "\b\n";
}

/// overload to show the content of this object
void holo_view_interactor::stream_stats(std::ostream& os)
{
	os << "holo_view_interactor:\n\a";

	oprintf(os, "y_view_angle=%.1fº, y_extent=%.1f, inp_z_range:[%.2f,%.2f]", y_view_angle, y_extent_at_focus, z_near,
			z_far);
	if (scene_extent.is_valid()) {
		oprintf(os, " adapted to scene: [%.2f,%.2f]\n", z_near_derived, z_far_derived);
		os << "current scene extent: " << scene_extent << std::endl;
	}
	else if (clip_relative_to_extent)
		oprintf(os, " adapted to extent: [%.2f,%.2f]\n", z_near_derived, z_far_derived);
	else
		os << "\n";

	oprintf(os, "foc=%.2f,%.2f,%.2f, dir=%.2f,%.2f,%.2f, up=%.2f,%.2f,%.2f\n", view::focus(0), view::focus(1),
			view::focus(2), view_dir(0), view_dir(1), view_dir(2), view_up_dir(0), view_up_dir(1), view_up_dir(2));

	oprintf(os, "render mpx:%s, eye-dist=%.3f", find_enum_name(HOLO_ENUMS, multiview_mpx_mode).c_str(), eye_distance);
	if (gamepad_attached) {
		oprintf(os, "\nleft:%s, right:%s", (left_mode == 0 ? "rotate" : "pan"),
				(right_mode == 0 ? "roll|zoom" : "move|dolly"));
	}
	if (emulation_active) {
		os << "\nEMU:";
		for (unsigned i = 0; i < 6; ++i) {
			os << " " << emulation_axes[i];
		}
	}
	os << "\b\n";
}

unsigned holo_view_interactor::get_viewport_index(unsigned col_index, unsigned row_index) const
{
	unsigned view_index = row_index * nr_viewport_columns + col_index;
	if (view_index < 0)
		view_index = 0;
	else if (view_index >= views.size())
		view_index = (unsigned)(views.size() - 1);

	return view_index;
}

/// call this function before a drawing process to support viewport splitting inside the draw call via the
/// activate/deactivate functions
void holo_view_interactor::enable_viewport_splitting(unsigned nr_cols, unsigned nr_rows)
{
	do_viewport_splitting = true;
	nr_viewport_columns = nr_cols;
	nr_viewport_rows = nr_rows;
	ensure_viewport_view_number(nr_cols * nr_rows);
	post_redraw();
}

/// check whether viewport splitting is activated and optionally set the number of columns and rows if corresponding
/// pointers are passed
bool holo_view_interactor::is_viewport_splitting_enabled(unsigned* nr_cols_ptr, unsigned* nr_rows_ptr) const
{
	if (do_viewport_splitting) {
		if (nr_cols_ptr)
			*nr_cols_ptr = nr_viewport_columns;
		if (nr_rows_ptr)
			*nr_rows_ptr = nr_viewport_rows;
	}
	return do_viewport_splitting;
}

/// disable viewport splitting
void holo_view_interactor::disable_viewport_splitting() { do_viewport_splitting = false; }

holo_view_interactor::ivec4 holo_view_interactor::split_viewport(const ivec4 vp, int col_idx, int row_idx) const
{
	ivec4 new_vp;
	new_vp[2] = vp[2] / nr_viewport_columns;
	new_vp[0] = vp[0] + col_idx * new_vp[2];
	new_vp[3] = vp[3] / nr_viewport_rows;
	new_vp[1] = vp[1] + row_idx * new_vp[3];
	if (viewport_shrinkage > 0) {
		new_vp[0] += viewport_shrinkage;
		new_vp[2] -= 2 * viewport_shrinkage;
		new_vp[1] += viewport_shrinkage;
		new_vp[3] -= 2 * viewport_shrinkage;
	}
	return new_vp;
}

/// inside the drawing process activate the sub-viewport with the given column and row indices, always terminate an
/// activated viewport with deactivate_split_viewport
void holo_view_interactor::activate_split_viewport(cgv::render::context& ctx, unsigned col_index, unsigned row_index)
{
	if (!do_viewport_splitting)
		return;
	const ivec4& current_vp = ctx.get_window_transformation_array().front().viewport;
	ivec4 new_vp = split_viewport(current_vp, col_index, row_index);
	ctx.push_window_transformation_array();
	ctx.set_viewport(new_vp);
	double aspect = (double)new_vp[2] / new_vp[3];
	unsigned view_index = get_viewport_index(col_index, row_index);
	ensure_viewport_view_number(view_index + 1);
	if (use_individual_view[view_index])
		compute_clipping_planes(views[view_index], z_near_derived, z_far_derived, clip_relative_to_extent);
	gl_set_projection_matrix(ctx, current_e, aspect);
	if (use_individual_view[view_index]) {
		compute_clipping_planes(z_near_derived, z_far_derived, clip_relative_to_extent);
		gl_set_modelview_matrix(ctx, current_e, aspect, views[view_index]);
	}
	((current_e == GLSU_RIGHT) ? MPWs_right : MPWs)[view_index] = ctx.get_modelview_projection_window_matrix();
}

/// deactivate the previously split viewport
void holo_view_interactor::deactivate_split_viewport(cgv::render::context& ctx)
{
	if (!do_viewport_splitting)
		return;
	ctx.pop_window_transformation_array();
	const ivec4& current_vp = ctx.get_window_transformation_array().front().viewport;
	double aspect = (double)current_vp[2] / current_vp[3];
	gl_set_projection_matrix(ctx, current_e, aspect);
	gl_set_modelview_matrix(ctx, current_e, aspect, *this);
}

/// make a viewport manage its own view
void holo_view_interactor::enable_viewport_individual_view(unsigned col_index, unsigned row_index, bool enable)
{
	unsigned view_index = get_viewport_index(col_index, row_index);
	ensure_viewport_view_number(view_index + 1);
	use_individual_view[view_index] = enable;
}

/// access the view of a given viewport
cgv::render::view& holo_view_interactor::ref_viewport_view(unsigned col_index, unsigned row_index)
{
	unsigned view_index = get_viewport_index(col_index, row_index);
	ensure_viewport_view_number(view_index + 1);
	return views[view_index];
}

//! given a mouse location and the pixel extent of the context, return the MPW matrix for unprojection
int holo_view_interactor::get_modelview_projection_window_matrices(
	  int x, int y, int width, int height, const dmat4** MPW_pptr, const dmat4** MPW_other_pptr, int* x_other_ptr,
	  int* y_other_ptr, int* vp_col_idx_ptr, int* vp_row_idx_ptr, int* vp_width_ptr, int* vp_height_ptr,
	  int* vp_center_x_ptr, int* vp_center_y_ptr, int* vp_center_x_other_ptr, int* vp_center_y_other_ptr) const
{
	*MPW_pptr = &MPW;
	const dmat4* MPW_other_ptr_local = &MPW;
	int vp_width = width;
	int vp_height = height;
	int eye_panel = 0;
	int off_x = 0, off_y = 0; // start of stereo viewport in mouse integer coordinates
	int off_x_other = 0, off_y_other = 0;
	int x_other = x;
	int y_other = y;
	int vp_col_idx = 0;
	int vp_row_idx = 0;
	if (last_do_viewport_splitting) {
		vp_width /= last_nr_viewport_columns;
		vp_height /= last_nr_viewport_rows;
		vp_col_idx = (x - off_x) / vp_width;
		vp_row_idx = (y - off_y) / vp_height;
		off_x += vp_col_idx * vp_width + viewport_shrinkage;
		off_y += vp_row_idx * vp_height + viewport_shrinkage;
		off_x_other += vp_col_idx * vp_width + viewport_shrinkage;
		off_y_other += vp_row_idx * vp_height + viewport_shrinkage;
		vp_width -= 2 * viewport_shrinkage;
		vp_height -= 2 * viewport_shrinkage;
		int vp_idx = vp_row_idx * last_nr_viewport_columns + vp_col_idx;
		if (eye_panel == 1) {
			if (vp_idx < (int)MPWs_right.size())
				*MPW_pptr = &MPWs_right[vp_idx];
			if (vp_idx < (int)MPWs.size())
				MPW_other_ptr_local = &MPWs[vp_idx];
		}
		else {
			if (vp_idx < (int)MPWs.size())
				*MPW_pptr = &MPWs[vp_idx];
			//			if (stereo_enabled) {
			//				if (vp_idx < (int)MPWs_right.size())
			//					MPW_other_ptr_local = &MPWs_right[vp_idx];
			//			}
			//			else {
			if (vp_idx < (int)MPWs.size())
				MPW_other_ptr_local = &MPWs[vp_idx];
			//			}
		}
	}

	if (MPW_other_pptr)
		*MPW_other_pptr = MPW_other_ptr_local;
	if (x_other_ptr)
		*x_other_ptr = x_other;
	if (y_other_ptr)
		*y_other_ptr = y_other;
	if (vp_col_idx_ptr)
		*vp_col_idx_ptr = vp_col_idx;
	if (vp_row_idx_ptr)
		*vp_row_idx_ptr = vp_row_idx;
	if (vp_width_ptr)
		*vp_width_ptr = vp_width;
	if (vp_height_ptr)
		*vp_height_ptr = vp_height;
	if (vp_center_x_ptr)
		*vp_center_x_ptr = off_x + vp_width / 2;
	if (vp_center_y_ptr)
		*vp_center_y_ptr = off_y + vp_height / 2;
	if (vp_center_x_other_ptr)
		*vp_center_x_other_ptr = off_x_other + vp_width / 2;
	if (vp_center_y_other_ptr)
		*vp_center_y_other_ptr = off_y_other + vp_height / 2;

	return eye_panel;
}

void holo_view_interactor::get_vp_col_and_row_indices(cgv::render::context& ctx, int x, int y, int& vp_col_idx,
													  int& vp_row_idx)
{
	const dmat4 *MPW_ptr, *MPW_other_ptr;
	int x_other, y_other, vp_width, vp_height;
	int eye_panel =
		  get_modelview_projection_window_matrices(x, y, ctx.get_width(), ctx.get_height(), &MPW_ptr, &MPW_other_ptr,
												   &x_other, &y_other, &vp_col_idx, &vp_row_idx, &vp_width, &vp_height);
}

double holo_view_interactor::get_z_and_unproject(cgv::render::context& ctx, int x, int y, dvec3& p)
{
	const dmat4 *MPW_ptr, *MPW_other_ptr;
	int x_other, y_other, vp_col_idx, vp_row_idx, vp_width, vp_height;
	int eye_panel =
		  get_modelview_projection_window_matrices(x, y, ctx.get_width(), ctx.get_height(), &MPW_ptr, &MPW_other_ptr,
												   &x_other, &y_other, &vp_col_idx, &vp_row_idx, &vp_width, &vp_height);
	ctx.make_current();
	double z = ctx.get_window_z(x, y);
	double z_other = ctx.get_window_z(x_other, y_other);

	if (z <= z_other) {
		p = ctx.get_model_point(x, y, z, *MPW_ptr);
		return z;
	}
	else {
		p = ctx.get_model_point(x_other, y_other, z_other, *MPW_other_ptr);
		return z_other;
	}
}

cgv::render::view::dvec3 unpack_dir(char c)
{
	switch (c) {
	case 'x':
		return cgv::render::view::dvec3(1, 0, 0);
	case 'X':
		return cgv::render::view::dvec3(-1, 0, 0);
	case 'y':
		return cgv::render::view::dvec3(0, 1, 0);
	case 'Y':
		return cgv::render::view::dvec3(0, -1, 0);
	case 'z':
		return cgv::render::view::dvec3(0, 0, 1);
	case 'Z':
		return cgv::render::view::dvec3(0, 0, -1);
	}
	return cgv::render::view::dvec3(0, 0, 0);
}

void holo_view_interactor::set_view_orientation(const std::string& axes)
{
	dvec3 axis;
	double angle;
	compute_axis_and_angle(unpack_dir(axes[0]), unpack_dir(axes[1]), axis, angle);
	cgv::gui::animate_with_axis_rotation(view_dir, axis, angle, 0.5)->set_base_ptr(this);
	cgv::gui::animate_with_axis_rotation(view_up_dir, axis, angle, 0.5)->set_base_ptr(this);
}
/// overload and implement this method to handle events
bool holo_view_interactor::handle(event& e)
{
	if (use_gamepad && ((e.get_flags() & EF_PAD) != 0)) {
		if (!gamepad_attached) {
			gamepad_attached = true;
			update_member(&gamepad_attached);
			post_redraw();
		}
		if (e.get_kind() == EID_THROTTLE) {
			cgv::gui::throttle_event& te = static_cast<cgv::gui::throttle_event&>(e);
			trigger[te.get_throttle_index()] = te.get_value();
			return true;
		}
		else if (e.get_kind() == EID_STICK) {
			cgv::gui::stick_event& se = static_cast<cgv::gui::stick_event&>(e);
			if (se.get_stick_index() == 0)
				left_stick = se.get_position();
			else
				right_stick = se.get_position();
			return true;
		}
	}
	if (e.get_kind() == EID_KEY) {
		key_event ke = (key_event&)e;
		if (gamepad_emulation) {
			switch (ke.get_key()) {
			case 'A':
				plus_key_action(0, ke.get_action());
				return true;
			case 'D':
				minus_key_action(0, ke.get_action());
				return true;
			case 'W':
				plus_key_action(1, ke.get_action());
				return true;
			case 'X':
				minus_key_action(1, ke.get_action());
				return true;
			case 'H':
				plus_key_action(2, ke.get_action());
				return true;
			case 'K':
				minus_key_action(2, ke.get_action());
				return true;
			case 'U':
				plus_key_action(3, ke.get_action());
				return true;
			case 'M':
				minus_key_action(3, ke.get_action());
				return true;
			}
		}
		if (ke.get_action() != KA_RELEASE) {
			switch (ke.get_key()) {
			case gamepad::GPK_LEFT_STICK_PRESS:
				left_mode = 1 - left_mode;
				on_set(&left_mode);
				return true;
			case gamepad::GPK_RIGHT_STICK_PRESS:
				right_mode = 1 - right_mode;
				on_set(&right_mode);
				return true;
			case gamepad::GPK_START:
				set_default_view();
				post_redraw();
				return true;
			case KEY_Space:
				if (ke.get_modifiers() == EM_CTRL) {
					set_default_view();
					post_redraw();
					return true;
				}
				break;
			case 'F':
			case gamepad::GPK_A:
				if (ke.get_modifiers() == EM_SHIFT || e.get_kind() == gamepad::GPK_A) {
					show_focus = !show_focus;
					on_set(&show_focus);
					return true;
				}
				else if (ke.get_modifiers() == EM_CTRL)
					z_far /= 1.05;
				else if (ke.get_modifiers() == EM_ALT)
					z_far *= 1.05;
				else
					break;
				on_set(&z_far);
				return true;
			case 'N':
				if (ke.get_modifiers() == EM_CTRL)
					z_near /= 1.05;
				else if (ke.get_modifiers() == EM_ALT)
					z_near *= 1.05;
				else
					break;
				on_set(&z_near);
				return true;
			case KEY_F4:
				if (ke.get_modifiers() == 0) {
					multiview_mpx_mode = MultiViewMode(((int&)multiview_mpx_mode + 1) % 3);
					on_set(&multiview_mpx_mode);
				}
				else
					break;
				return true;
			case KEY_Num_Add:
			case KEY_Page_Up:
				if ((ke.get_key() == KEY_Page_Up && ke.get_modifiers() == 0) ||
					(ke.get_key() == KEY_Num_Add && ke.get_modifiers() == EM_CTRL))
				{
					y_extent_at_focus /= pow(1.2, 1 / zoom_sensitivity);
					on_set(&y_extent_at_focus);
					return true;
				}
				break;
			case KEY_Num_Sub:
			case KEY_Page_Down:
				if ((ke.get_key() == KEY_Page_Down && ke.get_modifiers() == 0) ||
					(ke.get_key() == KEY_Num_Sub && ke.get_modifiers() == EM_CTRL))
				{
					y_extent_at_focus *= pow(1.2, 1 / zoom_sensitivity);
					on_set(&y_extent_at_focus);
					return true;
				}
				break;
			case gamepad::GPK_DPAD_RIGHT:
				set_view_orientation("xy");
				return true;
			case gamepad::GPK_DPAD_LEFT:
				set_view_orientation("Xy");
				return true;
			case gamepad::GPK_DPAD_UP:
				set_view_orientation("yz");
				return true;
			case gamepad::GPK_DPAD_DOWN:
				set_view_orientation("Yz");
				return true;
			case gamepad::GPK_LEFT_BUMPER:
				set_view_orientation("zy");
				return true;
			case gamepad::GPK_RIGHT_BUMPER:
				set_view_orientation("Zy");
				return true;
			case 'I':
				if (ke.get_modifiers() == cgv::gui::EM_SHIFT + cgv::gui::EM_CTRL) {
					quilt_interpolate = !quilt_interpolate;
					on_set(&quilt_interpolate);
					return true;
				}
				break;
			case 'S':
				if (ke.get_modifiers() == cgv::gui::EM_SHIFT + cgv::gui::EM_CTRL) {
					holo_storage_mode = HM_SINGLE;
					on_set(&holo_storage_mode);
					return true;
				}
				break;
			case 'V':
				if (ke.get_modifiers() == cgv::gui::EM_SHIFT + cgv::gui::EM_CTRL) {
					holo_storage_mode = HM_VOLUME;
					on_set(&holo_storage_mode);
					return true;
				}
				break;
			case 'Q':
				if (ke.get_modifiers() == cgv::gui::EM_SHIFT + cgv::gui::EM_CTRL) {
					holo_storage_mode = HM_QUILT;
					on_set(&holo_storage_mode);
					return true;
				}
				break;
			case 'X':
				if (ke.get_modifiers() == (cgv::gui::EM_SHIFT | cgv::gui::EM_CTRL))
					set_view_orientation("Xy");
				else if (ke.get_modifiers() == cgv::gui::EM_CTRL)
					set_view_orientation("xy");
				else
					break;
				return true;
			case 'Y':
				if (ke.get_modifiers() == (cgv::gui::EM_SHIFT | cgv::gui::EM_CTRL))
					set_view_orientation("Yz");
				else if (ke.get_modifiers() == cgv::gui::EM_CTRL)
					set_view_orientation("yz");
				else
					break;
				return true;
			case 'Z':
				if (ke.get_modifiers() == (cgv::gui::EM_SHIFT | cgv::gui::EM_CTRL))
					set_view_orientation("Zy");
				else if (ke.get_modifiers() == cgv::gui::EM_CTRL)
					set_view_orientation("zy");
				else
					break;
				return true;
			}
		}
	}
	else if (e.get_kind() == EID_MOUSE) {
		cgv::gui::mouse_event me = (cgv::gui::mouse_event&)e;
		int x_gl = me.get_x();
		int y_gl = get_context()->get_height() - 1 - me.get_y();
		if (me.get_action() == cgv::gui::MA_LEAVE)
			last_x = -1;
		else {
			last_x = x_gl;
			last_y = y_gl;
		}
		int width = 640, height = 480;
		int center_x = 320, center_y = 240;
		int vp_col_idx, vp_row_idx;
		cgv::render::view* view_ptr = this;
		const dmat4* MPW_ptr = 0;
		if (get_context()) {
			int eye = get_modelview_projection_window_matrices(
				  x_gl, y_gl, get_context()->get_width(), get_context()->get_height(), &MPW_ptr, 0, 0, 0, &vp_col_idx,
				  &vp_row_idx, &width, &height, &center_x, &center_y);
			unsigned view_index = get_viewport_index(vp_col_idx, vp_row_idx);
			if (view_index != -1 && use_individual_view[view_index]) {
				view_ptr = &views[view_index];
			}
		}
		dvec3 x, y, z;
		view_ptr->put_coordinate_system(x, y, z);

		switch (me.get_action()) {
		case MA_PRESS:
			if (me.get_button() == MB_LEFT_BUTTON && me.get_modifiers() == 0) {
				check_for_click = me.get_time();
				return true;
			}
			if (((me.get_button() == MB_LEFT_BUTTON) &&
				 ((me.get_modifiers() == 0) || (me.get_modifiers() == EM_SHIFT))) ||
				((me.get_button() == MB_RIGHT_BUTTON) && (me.get_modifiers() == 0)) ||
				((me.get_button() == MB_MIDDLE_BUTTON) && (me.get_modifiers() == 0)))
				return true;
			break;
		case MA_RELEASE:
			if (check_for_click != -1) {
				double dt = me.get_time() - check_for_click;
				if (dt < 0.2) {
					if (get_context()) {
						cgv::render::context& ctx = *get_context();
						dvec3 p;
						double z = get_z_and_unproject(ctx, x_gl, y_gl, p);
						if (z > 0 && z < 1) {
							if (y_view_angle > 0.1) {
								dvec3 e = view_ptr->get_eye();
								double l_old = (e - view_ptr->get_focus()).length();
								double l_new = dot(p - e, view_ptr->get_view_dir());
								// std::cout << "e=(" << e << "), p=(" << p << "), vd=(" << view_ptr->get_view_dir() <<
								// ") l_old=" << l_old << ", l_new=" << l_new << std::endl;
								cgv::gui::animate_with_geometric_blend(
									  view_ptr->ref_y_extent_at_focus(),
									  view_ptr->get_y_extent_at_focus() * l_new / l_old, 0.5)
									  ->set_base_ptr(this);
							}
							cgv::gui::animate_with_linear_blend(view_ptr->ref_focus(), p, 0.5)
								  ->configure(cgv::gui::APM_SIN_SQUARED, this);

							update_vec_member(view::focus);
							post_redraw();
							return true;
						}
					}
				}
				check_for_click = -1;
			}
			if ((me.get_button() == MB_LEFT_BUTTON && (me.get_modifiers() == 0 || me.get_modifiers() == EM_SHIFT)) ||
				me.get_button() == MB_RIGHT_BUTTON && me.get_modifiers() == 0)
				return true;
			break;
		case MA_MOVE:
			break;
		case MA_DRAG:
			check_for_click = -1;
			if (me.get_dx() == 0 && me.get_dy() == 0)
				break;
			if (me.get_button_state() == MB_LEFT_BUTTON && me.get_modifiers() == 0) {
				if (!two_d_enabled) {
					view_ptr->rotate(-6.0 * me.get_dy() / height / rotate_sensitivity,
									 -6.0 * me.get_dx() / width / rotate_sensitivity, view_ptr->get_depth_of_focus());
					update_vec_member(view_up_dir);
					update_vec_member(view_dir);
					post_redraw();
					return true;
				}
			}
			if (me.get_button_state() == MB_LEFT_BUTTON && me.get_modifiers() == EM_SHIFT) {
				int rx = me.get_x() - center_x;
				int ry = me.get_y() - center_y;
				double ds =
					  sqrt(((double)me.get_dx() * (double)me.get_dx() + (double)me.get_dy() * (double)me.get_dy()) /
						   ((double)rx * (double)rx + (double)ry * (double)ry));
				if (rx * me.get_dy() > ry * me.get_dx())
					ds = -ds;
				view_ptr->roll(ds / rotate_sensitivity);
				update_vec_member(view_up_dir);
				post_redraw();
				return true;
			}
			if (me.get_button_state() == MB_RIGHT_BUTTON && me.get_modifiers() == 0) {
				view_ptr->set_focus(view_ptr->get_focus() -
									(view_ptr->get_y_extent_at_focus() * me.get_dx() / width) * x +
									(view_ptr->get_y_extent_at_focus() * me.get_dy() / height) * y);
				update_vec_member(view::focus);
				post_redraw();
				return true;
			}
			if (me.get_button_state() == MB_MIDDLE_BUTTON && me.get_modifiers() == 0) {
				view_ptr->set_focus(view_ptr->get_focus() - 5 * view_ptr->get_y_extent_at_focus() * me.get_dy() /
																  height * z / zoom_sensitivity);
				update_vec_member(view::focus);
				post_redraw();
				return true;
			}
			break;
		case MA_WHEEL:
			if (e.get_modifiers() == EM_ALT) {
				eye_distance -= 0.001 * me.get_dy();
				if (eye_distance < 0)
					eye_distance = 0;
				update_member(&eye_distance);
				post_redraw();
				return true;
			}
			else if (e.get_modifiers() == EM_CTRL) {
				parallax_zero_scale *= exp(-0.03 * me.get_dy());
				if (parallax_zero_scale > 1)
					parallax_zero_scale = 1;
				else if (parallax_zero_scale < 0.01)
					parallax_zero_scale = 0.01;
				update_member(&parallax_zero_scale);
				post_redraw();
				return true;
			}
			else if (e.get_modifiers() == EM_SHIFT) {
				view_ptr->set_y_view_angle(view_ptr->get_y_view_angle() + me.get_dy() * 5);
				if (view_ptr->get_y_view_angle() < 0)
					view_ptr->set_y_view_angle(0);
				if (view_ptr->get_y_view_angle() > 180)
					view_ptr->set_y_view_angle(180);
				update_member(&y_view_angle);
				post_redraw();
				return true;
			}
			else if (e.get_modifiers() == 0) {
				double scale = exp(0.2 * me.get_dy() / zoom_sensitivity);
				view_ptr->set_y_extent_at_focus(view_ptr->get_y_extent_at_focus() * scale);
				update_member(&y_extent_at_focus);
				post_redraw();
				return true;
			}
			break;
		default:
			break;
		}
	}
	return false;
}

void holo_view_interactor::on_rotation_change()
{
	for (unsigned i = 0; i < 3; ++i) {
		if (fix_view_up_dir)
			view_up_dir(i) = i == 1 ? 1 : 0;
		update_member(&view_up_dir(i));
		update_member(&view_dir(i));
	}
	post_redraw();
}
/// this method is called in one pass over all drawables after drawing
void holo_view_interactor::finish_frame(cgv::render::context& ctx)
{
	cgv::render::RenderPassFlags rpf = ctx.get_render_pass_flags();
	if ((rpf & RPF_SET_MODELVIEW_PROJECTION) == 0)
		return;
	if (show_focus) {
		/*		ctx.push_P();
				ctx.push_V();
				ctx.set_P(P);
				ctx.set_V(V);
				*/
		glDisable(GL_DEPTH_TEST);
		glLineStipple(1, 15);
		glEnable(GL_LINE_STIPPLE);

		if (is_viewport_splitting_enabled()) {
			for (unsigned c = 0; c < nr_viewport_columns; ++c) {
				for (unsigned r = 0; r < nr_viewport_rows; ++r) {
					activate_split_viewport(ctx, c, r);
					draw_focus();
					deactivate_split_viewport(ctx);
				}
			}
		}
		else
			draw_focus();
		glDisable(GL_LINE_STIPPLE);
		glEnable(GL_DEPTH_TEST);
		/*
		ctx.pop_V();
		ctx.pop_P();
		*/
	}
}

/// set the current projection matrix
void holo_view_interactor::gl_set_projection_matrix(cgv::render::context& ctx, float e, double aspect)
{
	dmat4 P;
	if (y_view_angle <= 0.1)
		P = ortho4<double>(-aspect * y_extent_at_focus, aspect * y_extent_at_focus, -y_extent_at_focus,
						   y_extent_at_focus, z_near_derived, z_far_derived);
	else {
		if (stereo_translate_in_model_view)
			P = cgv::math::stereo_frustum_screen4<double>(e, eye_distance, y_extent_at_focus * aspect,
														  y_extent_at_focus, get_parallax_zero_depth(), z_near_derived,
														  z_far_derived);
		else
			P = cgv::math::stereo_perspective_screen4<double>(e, eye_distance, y_extent_at_focus * aspect,
															  y_extent_at_focus, get_parallax_zero_depth(),
															  z_near_derived, z_far_derived);
	}
	ctx.set_projection_matrix(P);
}

void holo_view_interactor::gl_set_modelview_matrix(cgv::render::context& ctx, float e, double aspect,
												   const cgv::render::view& view)
{
	ctx.set_modelview_matrix(cgv::math::identity4<double>());
	if (stereo_translate_in_model_view)
		ctx.mul_modelview_matrix(
			  cgv::math::stereo_translate_screen4<double>(e, eye_distance, view.get_y_extent_at_focus() * aspect));
	ctx.mul_modelview_matrix(cgv::math::look_at4(view.get_eye(), view.get_focus(), view.get_view_up_dir()));
}

/// ensure sufficient number of viewport views
void holo_view_interactor::ensure_viewport_view_number(unsigned nr)
{
	if (views.size() < nr) {
		unsigned old_nr = (unsigned)views.size();
		views.resize(nr);
		use_individual_view.resize(nr);
		for (unsigned i = old_nr; i < nr; ++i) {
			views[i] = *this;
			use_individual_view[i] = false;
		}
	}
}

bool holo_view_interactor::init(cgv::render::context& ctx)
{
	// create offscreen framebuffers used for warping approaches
	for (unsigned i = 0; i < 3; i++) {
		render_fbo[i].add_attachment("depth", "[D]");
		render_fbo[i].add_attachment("color", "uint8[R,G,B,A]");
	}

	// add shader programs to a shader library and load all of them
	if (!shader_lib.add("quilt_prog", "quilt_finish.glpr"))
		return false;
	if (!shader_lib.add("volume_prog", "volume_finish.glpr"))
		return false;
	if (!shader_lib.add("reproject_shader", "reprojection.glpr"))
		return false;
	if (!shader_lib.add("vwarp_shader", "vertex_warp.glpr"))
		return false;
	if (!shader_lib.add("compute_shader", "compute_warp.glpr"))
		return false;
	if (!shader_lib.add("volume_resolve_compute_shader", "compute_resolve_volume.glpr"))
		return false;
	if (!shader_lib.add("quilt_resolve_compute_shader", "compute_resolve_quilt.glpr"))
		return false;
	if (!shader_lib.load_all(ctx))
		return false;

	// set shader programs
	quilt_prog = shader_lib.get("quilt_prog");
	volume_prog = shader_lib.get("volume_prog");
	reproject_shader = shader_lib.get("reproject_shader");
	vwarp_shader = shader_lib.get("vwarp_shader");
	compute_shader = shader_lib.get("compute_shader");

	view_width = ctx.get_width();
	view_height = ctx.get_height();

	//generate shader storage buffer for compute shader warping
	glGenBuffers(1, &ssbo);
	glBindBuffer(GL_ARRAY_BUFFER, ssbo);

	// set up size of ssbo 
	glNamedBufferData(ssbo,
					  GLsizeiptr(sizeof(unsigned int) * view_width * view_height *
								 std::ceil(nr_holo_views / (float)quilt_nr_cols) * quilt_nr_cols),
						nullptr, GL_DYNAMIC_COPY);

	#if EVAL == 1
	// generate time query for evaluation
	glGenQueries(1, &time_query);

	// in case the .csv file for the performance measurements is empty or doesn't exist, set it up
	std::ifstream in("measurements.csv");
	if (in.is_open()) {
		in.seekg(0, std::ios::end);
		size_t size = in.tellg();
		if (size == 0) {
			set_up_eval_file();
		}
	}
	else {
		set_up_eval_file();
	}
	#endif

	return true;
}

/// this method is called in one pass over all drawables before the draw method
void holo_view_interactor::init_frame(context& ctx)
{
	cgv::render::RenderPassFlags rpf = ctx.get_render_pass_flags();

	//set size of quilt so that all views fit into it
	quilt_width = view_width * quilt_nr_cols;
	quilt_height = view_height * quilt_nr_rows;

	if (!mesh)
		mesh = &dynamic_cast<mesh_viewer&>(*cgv::base::find_object_by_name("mesh_viewer"));

	// update resolve compute shaders depending on whether or not splatting should be done
	if (update_defines) {
		if (multiview_mpx_mode == MVM_COMPUTE)
			compute_define["SPLAT"] = "0";
		else
			compute_define["SPLAT"] = "1";
		shader_lib.reload(ctx, "volume_resolve_compute_shader", compute_define);
		shader_lib.reload(ctx, "quilt_resolve_compute_shader", compute_define);
		update_defines = false;
	}

	switch (multiview_mpx_mode) {
	case MVM_SINGLE:
		current_e = (2.0f * view_index) / (nr_holo_views - 1) - 1.0f;
		last_do_viewport_splitting = do_viewport_splitting;
		last_nr_viewport_columns = nr_viewport_columns;
		last_nr_viewport_rows = nr_viewport_rows;
		break;
	case MVM_CONVENTIONAL:
		// render the amount of holo views in individual render passes
		if (initiate_render_pass_recursion(ctx)) {
			enable_surface(ctx);
			last_do_viewport_splitting = do_viewport_splitting;
			last_nr_viewport_columns = nr_viewport_columns;
			last_nr_viewport_rows = nr_viewport_rows;
			if (holo_storage_mode == HM_QUILT) {
				ctx.set_bg_color(quilt_bg_color.R(), quilt_bg_color.G(), quilt_bg_color.B(), 1.0f);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				vi = 0;
				for (quilt_row = 0; quilt_row < quilt_nr_rows; ++quilt_row) {
					for (quilt_col = 0; quilt_col < quilt_nr_cols; ++quilt_col) {
						volume_fbo.attach(ctx, volume_holo_tex, view_index, 0, 0);
						perform_render_pass(ctx, vi, RP_STEREO);
						if (++vi == nr_holo_views)
							break;
					}
					if (vi == nr_holo_views)
						break;
				}
			}
			else {
				ctx.set_bg_color(0, 0, 0, 1);
				for (vi = 0; vi < nr_holo_views; ++vi) {
					volume_fbo.attach(ctx, volume_holo_tex, vi, 0, 0);
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
					perform_render_pass(ctx, vi, RP_STEREO);
				}
			}
			initiate_terminal_render_pass(nr_holo_views - 1);
		}
		if (!multi_pass_ignore_finish(ctx)) {
			current_e = (2.0f * vi) / (nr_holo_views - 1) - 1.0f;
			if (holo_storage_mode == HM_QUILT) {
				ivec4 vp(quilt_col * view_width, quilt_row * view_height, view_width, view_height);
				glViewport(vp[0], vp[1], vp[2], vp[3]);
				glScissor(vp[0], vp[1], vp[2], vp[3]);
				glEnable(GL_SCISSOR_TEST);
			}
		}
		break;
	case MVM_GEOMETRY:
		// render the amount of holo views in render passes using a geometry shader where four views are always rendered at once
		if (initiate_render_pass_recursion(ctx)) {
			enable_surface(ctx);
			last_do_viewport_splitting = do_viewport_splitting;
			last_nr_viewport_columns = nr_viewport_columns;
			last_nr_viewport_rows = nr_viewport_rows;
			if (holo_storage_mode == HM_QUILT) {
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				ctx.set_bg_color(quilt_bg_color.R(), quilt_bg_color.G(), quilt_bg_color.B(), 1.0f);
				vi = 0;
				while (vi < nr_holo_views) {
					volume_fbo.attach(ctx, volume_holo_tex, view_index, 0, 0);
					// already compute current_e as it is needed for the geometry pass
					current_e = (2.0f * vi) / (nr_holo_views - 1) - 1.0f;
					mesh->set_params_for_gemoetry(get_parallax_zero_depth(), eye_distance, current_e,
												  (float)nr_holo_views);
					perform_render_pass(ctx, vi, RP_STEREO);
				}
			}
			else {
				vi = 0;
				while (vi < nr_holo_views) {
					// attach layered framebuffer and textures 
					glBindFramebuffer(GL_FRAMEBUFFER, (unsigned)((size_t)layered_fbo.handle) - 1);
					glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
										 (unsigned)((size_t)layered_depth_tex.handle) - 1, 0);
					glClear(GL_COLOR_BUFFER_BIT);
					glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
										 (unsigned)((size_t)layered_color_tex.handle) - 1, 0);
					glClear(GL_COLOR_BUFFER_BIT);

					// already compute current_e as it is needed for the geometry pass
					current_e = (2.0f * vi) / (nr_holo_views - 1) - 1.0f;
					mesh->set_params_for_gemoetry(get_parallax_zero_depth(), eye_distance, current_e,
												  (float)nr_holo_views);
					perform_render_pass(ctx, vi, RP_STEREO);

					// copy result from layered texture to volume texture
					for (int i = 0; i < 4; i++) {
						glCopyImageSubData((unsigned)((size_t)layered_color_tex.handle) - 1, GL_TEXTURE_2D_ARRAY, 0, 0,
										   0, i, (unsigned)((size_t)volume_holo_tex.handle) - 1, GL_TEXTURE_3D, 0, 0, 0,
										   vi, view_width, view_height, 1);
						if (++vi == nr_holo_views)
							break;
					}
					if (vi == nr_holo_views)
						break;

					glBindFramebuffer(GL_FRAMEBUFFER, 0);
				}
			}
			initiate_terminal_render_pass(nr_holo_views - 1);
		}
		if (!multi_pass_ignore_finish(ctx)) {
			if (holo_storage_mode == HM_QUILT) {
				glEnable(GL_SCISSOR_TEST);
				// set the next four viewports for the next geoemtry shader render pass in case of quilt mode
				for (GLuint i = 0; i < 4; i++) {
					quilt_row = vi / quilt_nr_cols;
					quilt_col = vi % quilt_nr_cols;
					ivec4 vp(quilt_col * view_width, quilt_row * view_height, view_width, view_height);
					glViewportIndexedf(i, vp[0], vp[1], vp[2], vp[3]);
					glScissorIndexed(i, vp[0], vp[1], vp[2], vp[3]);
					vi++;
				}
			}
		}
		break;
	case MVM_REPROJECT:
	case MVM_VWARP:
	case MVM_VWARP_CLOSEST:
	case MVM_COMPUTE:
	case MVM_COMPUTE_SPLAT:
		// render one, two, or three source views into offscreen framebuffer for all warping approaches
		if (initiate_render_pass_recursion(ctx)) {
			enable_surface(ctx);
			last_do_viewport_splitting = do_viewport_splitting;
			last_nr_viewport_columns = nr_viewport_columns;
			last_nr_viewport_rows = nr_viewport_rows;
			vi = 0;
			for (int i = 0; i < nr_render_views; i++) {
				current_render_fbo = render_fbo[vi];
				current_render_fbo.enable(ctx);
				current_render_fbo.ref_frame_buffer().push_viewport(ctx);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

				perform_render_pass(ctx, vi, RP_STEREO);

				vi++;
			}
			initiate_terminal_render_pass(nr_render_views - 1);
		}
		if (!multi_pass_ignore_finish(ctx)) {
			current_e = nr_render_views == 1 ? 0.0 : (2.0f * vi) / (nr_render_views - 1) - 1.0f;
		}
		break;
	}

	double aspect = (double)view_width / view_height;

	// compute the clipping planes based on the eye and scene extent
	compute_clipping_planes(z_near_derived, z_far_derived, clip_relative_to_extent);
	if (rpf & RPF_SET_PROJECTION) {
		gl_set_projection_matrix(ctx, current_e, aspect);
		// store projection matrix for warping approaches and calculate position of source views
		if ((vi < nr_render_views) && multiview_mpx_mode != MVM_CONVENTIONAL && multiview_mpx_mode != MVM_GEOMETRY) {
			proj_source[vi] = ctx.get_projection_matrix();
			eye_source[vi] = vec4(0.5f * current_e * eye_distance * y_extent_at_focus * aspect, 0, 0, 1);
		}
	}

	if (rpf & RPF_SET_MODELVIEW) {
		gl_set_modelview_matrix(ctx, current_e, aspect, *this);
		// store modelview matrix for reprojection 
		if ((vi < nr_render_views) && multiview_mpx_mode == MVM_REPROJECT)
			modelview_source[vi] = ctx.get_modelview_matrix();
	}

	if (current_e == GLSU_RIGHT) {
		MPW_right = ctx.get_modelview_projection_window_matrix();
		if (do_viewport_splitting)
			MPWs_right = std::vector<dmat4>(nr_viewport_rows * nr_viewport_columns, MPW_right);
	}
	else {
		MPW = ctx.get_modelview_projection_window_matrix();
		if (do_viewport_splitting)
			MPWs = std::vector<dmat4>(nr_viewport_rows * nr_viewport_columns, MPW);
	}
}

/// this method is called in one pass over all drawables after finish frame
void holo_view_interactor::after_finish(cgv::render::context& ctx)
{
	if (multiview_mpx_mode == MVM_CONVENTIONAL || multiview_mpx_mode == MVM_GEOMETRY) {
		if (!multi_pass_ignore_finish(ctx) && multi_pass_terminate(ctx)) {
			disable_surface(ctx);
			glScissor(0, 0, ctx.get_width(), ctx.get_height());
			glDisable(GL_SCISSOR_TEST);
			post_process_surface(ctx);
		}
	}
	else if (multiview_mpx_mode != MVM_SINGLE) {
		if (!multi_pass_ignore_finish(ctx) && multi_pass_terminate(ctx)) {
			// warp source views to target views
			compute_holo_views(ctx);
			disable_surface(ctx);
			post_process_surface(ctx);
		}
		// disable just used offscreen framebuffer
		else {
			current_render_fbo.ref_frame_buffer().pop_viewport(ctx);
			current_render_fbo.disable(ctx);
		}
	}
}

// create framebuffers, textures, and buffers, and enable the needed framebuffer
void holo_view_interactor::enable_surface(cgv::render::context& ctx)
{
	#if EVAL == 1
	// start evaluation 
	if (evaluate)
		glBeginQuery(GL_TIME_ELAPSED, time_query);
	#endif

	if (holo_storage_mode == HM_QUILT) {
		if (!quilt_use_offline_texture)
			return;
		if (!quilt_fbo.is_created() || quilt_fbo.get_width() != quilt_width ||
			quilt_fbo.get_height() != quilt_height)
		{
			quilt_holo_tex.destruct(ctx);
			quilt_fbo.destruct(ctx);
			quilt_depth_buffer.destruct(ctx);
			quilt_fbo.create(ctx, quilt_width, quilt_height);
			quilt_depth_buffer.create(ctx, quilt_width, quilt_height);
			quilt_holo_tex.create(ctx, TT_2D, quilt_width, quilt_height);
			quilt_fbo.attach(ctx, quilt_holo_tex);
			quilt_fbo.attach(ctx, quilt_depth_buffer);
		}
		quilt_fbo.enable(ctx);
		quilt_fbo.push_viewport(ctx);
	}
	else {
		if (!volume_fbo.is_created() || volume_fbo.get_width() != view_width ||
			volume_fbo.get_height() != view_height || change_nr_holo_views)
		{
			volume_fbo.destruct(ctx);
			volume_holo_tex.destruct(ctx);
			volume_depth_buffer.destruct(ctx);
			volume_fbo.create(ctx, view_width, view_height);
			volume_depth_buffer.create(ctx, view_width, view_height);
			volume_holo_tex.create(ctx, TT_3D, view_width, view_height, nr_holo_views);
			volume_fbo.attach(ctx, volume_holo_tex, 0, 0, 0);
			volume_fbo.attach(ctx, volume_depth_buffer);
			change_nr_holo_views = false;
		}
		volume_fbo.enable(ctx);
		volume_fbo.push_viewport(ctx);

		// create a layered framebuffer with depth of four in case the volume storage mode is used for the geometry approach
		if ((!layered_fbo.is_created() || layered_fbo.get_width() != view_width ||
			 layered_fbo.get_height() != view_height || layered_depth_tex.get_depth() != 4) &&
			multiview_mpx_mode == MVM_GEOMETRY && holo_storage_mode != HM_QUILT)
		{
			layered_fbo.destruct(ctx);
			layered_depth_tex.destruct(ctx);
			layered_color_tex.destruct(ctx);
			layered_depth_tex.create(ctx, TT_2D_ARRAY, view_width, view_height, 4);
			layered_color_tex.create(ctx, TT_2D_ARRAY, view_width, view_height, 4);
			layered_fbo.create(ctx, view_width, view_height);
		}
	}

	// when warping, ensure that offscreen framebuffers have the correct size
	if (multiview_mpx_mode == MVM_REPROJECT || multiview_mpx_mode == MVM_VWARP ||
		multiview_mpx_mode == MVM_VWARP_CLOSEST || multiview_mpx_mode == MVM_COMPUTE || multiview_mpx_mode == MVM_COMPUTE_SPLAT)
	{
		auto& fb = render_fbo[1].ref_frame_buffer();
		if (fb.get_width() != (int)view_width || fb.get_height() != (int)view_height) {
			render_fbo[0].set_size({(int)view_width, (int)view_height});
			render_fbo[0].ensure(ctx);

			const auto res = render_fbo[0].get_size();
			render_fbo[1].set_size(res);
			render_fbo[1].ensure(ctx);
			render_fbo[2].set_size(res);
			render_fbo[2].ensure(ctx);

			const float half_aspect = (float)res.x() / (2 * res.y());
			// generate heightmaps by tessellating the screen
			heightmap_reproject =
				  tessellator::quad(ctx, reproject_shader, {-half_aspect, -.5f, .0f},
												   {half_aspect, .5f, .0f}, res.x(), res.y(), tessellator::VA_TEXCOORD);
			heightmap_vwarp = tessellator::quad(ctx, vwarp_shader, {-half_aspect, -.5f, .0f}, {half_aspect, .5f, .0f},
												res.x(), res.y(), tessellator::VA_TEXCOORD);
		}
	}
}

void holo_view_interactor::disable_surface(cgv::render::context& ctx)
{
	if (holo_storage_mode == HM_QUILT) {
		quilt_fbo.pop_viewport(ctx);
		quilt_fbo.disable(ctx);
		// write the quilt to a file
		if (quilt_write_to_file) {
			std::filesystem::path cwd = std::filesystem::current_path();
			if (cwd.string().find("res") != std::string::npos)
				quilt_holo_tex.write_to_file(ctx, "../quilt.png");
			else
				quilt_holo_tex.write_to_file(ctx, "quilt.png");
			quilt_write_to_file = false;
			on_set(&quilt_write_to_file);
		}
	}
	else {
		// write current view index to a file
		if (volume_write_to_file) {
			std::filesystem::path cwd = std::filesystem::current_path();
			if (cwd.string().find("res") != std::string::npos)
				volume_holo_tex.write_to_file(ctx, "../volume.png", view_index);
			else
				volume_holo_tex.write_to_file(ctx, "volume.png", view_index);
			volume_write_to_file = false;
			on_set(&volume_write_to_file);
		}
		volume_fbo.pop_viewport(ctx);
		volume_fbo.disable(ctx);
	}
	glViewport(0, 0, ctx.get_width(), ctx.get_height());
}

// iterate once over all source views to generate one holo view with the reprojection approach
void holo_view_interactor::draw_reproject(cgv::render::context& ctx)
{
	reproject_shader.enable(ctx);
	glDisable(GL_CULL_FACE);
	reproject_shader.set_uniform(ctx, "epsilon", epsilon);
	reproject_shader.set_uniform(ctx, "artefacts", dis_artefacts);

	for (unsigned int i = 0; i < nr_render_views; i++) {
		texture &color_tex = *render_fbo[i].attachment_texture_ptr("color"),
				&depth_tex = *render_fbo[i].attachment_texture_ptr("depth");

		color_tex.enable(ctx, 0);
		reproject_shader.set_uniform(ctx, "color", 0);
		depth_tex.enable(ctx, 1);
		reproject_shader.set_uniform(ctx, "depth", 1);
		reproject_shader.set_uniform(ctx, "inv_mvp_source", inv(proj_source[i] * modelview_source[i]));

		heightmap_reproject.draw(ctx);
		color_tex.disable(ctx);
		depth_tex.disable(ctx);
	}
	glEnable(GL_CULL_FACE);
	reproject_shader.disable(ctx);
}

// generate one holo view with the vertex warping approach by warping the closest of the source views
void holo_view_interactor::draw_vertex_warp_closest(cgv::render::context& ctx)
{
	double aspect = (double)view_width / view_height;
	vec4 eye_target = vec4(0.5f * current_e * eye_distance * y_extent_at_focus * aspect, 0, 0, 1);

	// find closest source view
	int source_idx = 0;
	float diff = std::numeric_limits<float>::max();
	for (int i = 0; i < nr_render_views; i++) {
		if (abs(eye_target[0] - eye_source[i][0]) < diff) {
			source_idx = i;
			diff = abs(eye_target[0] - eye_source[i][0]);
		}
	}
	texture &color_tex = *render_fbo[source_idx].attachment_texture_ptr("color"),
			&depth_tex = *render_fbo[source_idx].attachment_texture_ptr("depth");

	float shear = (eye_target[0] - eye_source[source_idx][0]) * (get_parallax_zero_depth() - z_far_derived) /
				  get_parallax_zero_depth();

	color_tex.enable(ctx, 0);
	vwarp_shader.set_uniform(ctx, "color_tex", 0);
	depth_tex.enable(ctx, 1);
	vwarp_shader.set_uniform(ctx, "depth_tex", 1);
	vwarp_shader.set_uniform(ctx, "p_source", proj_source[source_idx]);
	vwarp_shader.set_uniform(ctx, "eye_source", eye_source[source_idx]);
	vwarp_shader.set_uniform(ctx, "eye_target", eye_target);
	vwarp_shader.set_uniform(ctx, "z_far", (float)z_far_derived);
	vwarp_shader.set_uniform(ctx, "shear", shear);
	vwarp_shader.set_uniform(ctx, "epsilon", epsilon);
	vwarp_shader.set_uniform(ctx, "artefacts", dis_artefacts);

	vwarp_shader.enable(ctx);
	glDisable(GL_CULL_FACE);
	heightmap_vwarp.draw(ctx);
	glEnable(GL_CULL_FACE);
	vwarp_shader.disable(ctx);
	color_tex.disable(ctx);
	depth_tex.disable(ctx);
}

// iterate once over all rendered views to generate one holo view with the vertex warping approach
void holo_view_interactor::draw_vertex_warp(cgv::render::context& ctx)
{
	double aspect = (double)view_width / view_height;
	vec4 eye_target = vec4(0.5f * current_e * eye_distance * y_extent_at_focus * aspect, 0, 0, 1);

	vwarp_shader.enable(ctx);
	glDisable(GL_CULL_FACE);

	vwarp_shader.set_uniform(ctx, "eye_target", eye_target);
	vwarp_shader.set_uniform(ctx, "z_far", (float)z_far_derived);
	vwarp_shader.set_uniform(ctx, "epsilon", epsilon);
	vwarp_shader.set_uniform(ctx, "artefacts", dis_artefacts);

	for (unsigned int i = 0; i < nr_render_views; i++) {
		texture &color_tex = *render_fbo[i].attachment_texture_ptr("color"),
				&depth_tex = *render_fbo[i].attachment_texture_ptr("depth");

		float shear = (eye_target[0] - eye_source[i][0]) * (get_parallax_zero_depth() - z_far_derived) /
					  get_parallax_zero_depth();

		color_tex.enable(ctx, 0);
		vwarp_shader.set_uniform(ctx, "color_tex", 0);
		depth_tex.enable(ctx, 1);
		vwarp_shader.set_uniform(ctx, "depth_tex", 1);
		vwarp_shader.set_uniform(ctx, "p_source", proj_source[i]);
		vwarp_shader.set_uniform(ctx, "eye_source", eye_source[i]);
		vwarp_shader.set_uniform(ctx, "shear", shear);

		heightmap_vwarp.draw(ctx);
		color_tex.disable(ctx);
		depth_tex.disable(ctx);
	}
	glEnable(GL_CULL_FACE);
	vwarp_shader.disable(ctx);
}

// do one compute shader call in order to warp all source views to all target views 
void holo_view_interactor::warp_compute_shader(cgv::render::context& ctx)
{
	double aspect = (double)view_width / view_height;
	float views_x_extent = eye_distance * y_extent_at_focus * aspect;
	float shear = (get_parallax_zero_depth() - z_far_derived) / get_parallax_zero_depth();

	texture &color_tex0 = *render_fbo[0].attachment_texture_ptr("color"),
			&depth_tex0 = *render_fbo[0].attachment_texture_ptr("depth"),
			&color_tex1 = *render_fbo[1].attachment_texture_ptr("color"),
			&depth_tex1 = *render_fbo[1].attachment_texture_ptr("depth"),
			&color_tex2 = *render_fbo[2].attachment_texture_ptr("color"),
			&depth_tex2 = *render_fbo[2].attachment_texture_ptr("depth");

	compute_shader.enable(ctx);

	// attach shader storage buffer object
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);

	// clear the ssbo with a value at the largest depth and an alpha value of zero (not visible)
	uint8_t clear_color[4] = {0, 0, 0, 254};
	glClearNamedBufferData(ssbo, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, &clear_color);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);

	color_tex0.enable(ctx, 0);
	compute_shader.set_uniform(ctx, "color_tex0", 0);
	depth_tex0.enable(ctx, 1);
	compute_shader.set_uniform(ctx, "depth_tex0", 1);
	color_tex1.enable(ctx, 2);
	compute_shader.set_uniform(ctx, "color_tex1", 2);
	depth_tex1.enable(ctx, 3);
	compute_shader.set_uniform(ctx, "depth_tex1", 3);
	color_tex2.enable(ctx, 4);
	compute_shader.set_uniform(ctx, "color_tex2", 4);
	depth_tex2.enable(ctx, 5);
	compute_shader.set_uniform(ctx, "depth_tex2", 5);

	compute_shader.set_uniform(ctx, "p_source_zero", proj_source[0]);
	compute_shader.set_uniform(ctx, "start_x", -views_x_extent / 2);
	compute_shader.set_uniform(ctx, "x_offset", views_x_extent / (nr_holo_views-1));
	compute_shader.set_uniform(ctx, "nr_holo_views", (int)nr_holo_views);
	compute_shader.set_uniform(ctx, "nr_render_views", (int)nr_render_views);
	compute_shader.set_uniform(ctx, "z_far", (float)z_far_derived);
	compute_shader.set_uniform(ctx, "shear", shear);
	compute_shader.set_uniform(ctx, "eye_sep", (float)eye_distance);
	compute_shader.set_uniform(ctx, "zero_parallax", (float)get_parallax_zero_depth());
	compute_shader.set_uniform(ctx, "screen_w", (int)view_width);
	compute_shader.set_uniform(ctx, "screen_h", (int)view_height);
	compute_shader.set_uniform(ctx, "quilt_cols", (int)quilt_nr_cols);

	// get work group size from compute shader
	glGetProgramiv((unsigned)((size_t)compute_shader.handle) - 1, GL_COMPUTE_WORK_GROUP_SIZE, local_work_group);
	// compute call 
	glDispatchCompute(ceil(view_width / (float)local_work_group[0]),
					  ceil(view_height / (float)local_work_group[1]), 1);

	glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);

	compute_shader.disable(ctx);
}

// compute shader to transfer data from shader storage buffer to volume texture
void holo_view_interactor::volume_resolve_pass_compute_shader(cgv::render::context& ctx)
{
	auto& volume_resolve_compute_shader = shader_lib.get("volume_resolve_compute_shader");
	volume_resolve_compute_shader.enable(ctx);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);

	glBindImageTexture(0, (int&)volume_holo_tex.handle - 1, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	volume_resolve_compute_shader.set_uniform(ctx, "screen_w", (int)view_width);
	volume_resolve_compute_shader.set_uniform(ctx, "screen_h", (int)view_height);
	volume_resolve_compute_shader.set_uniform(ctx, "quilt_cols", (int)quilt_nr_cols);
	volume_resolve_compute_shader.set_uniform(ctx, "nr_holo_views", (int)nr_holo_views);

	// get work group size from compute shader
	glGetProgramiv((unsigned)((size_t)volume_resolve_compute_shader.handle) - 1, GL_COMPUTE_WORK_GROUP_SIZE,
				   local_work_group);
	// compute call
	glDispatchCompute(ceil(view_width / (float)local_work_group[0]),
					  ceil(view_height / (float)local_work_group[1]), 1);

	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
	glBindImageTexture(0, 0, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	volume_resolve_compute_shader.disable(ctx);
}

// compute shader to transfer data from shader storage buffer to quilt texture
void holo_view_interactor::quilt_resolve_pass_compute_shader(cgv::render::context& ctx)
{
	auto& quilt_resolve_compute_shader = shader_lib.get("quilt_resolve_compute_shader");
	quilt_resolve_compute_shader.enable(ctx);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);

	glBindImageTexture(0, (int&)quilt_holo_tex.handle - 1, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	quilt_resolve_compute_shader.set_uniform(ctx, "screen_w", (int)view_width);
	quilt_resolve_compute_shader.set_uniform(ctx, "quilt_cols", (int)quilt_nr_cols);

	// get work group size from compute shader
	glGetProgramiv((unsigned)((size_t)quilt_resolve_compute_shader.handle) - 1, GL_COMPUTE_WORK_GROUP_SIZE,
				   local_work_group);
	// compute call
	glDispatchCompute(ceil(view_width * quilt_nr_cols / (float)local_work_group[0]),
					  ceil(view_height * quilt_nr_rows / (float)local_work_group[1]), 1);

	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
	glBindImageTexture(0, 0, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	quilt_resolve_compute_shader.disable(ctx);
}

// generate all needed target views out of the source views depending on what method is used
void holo_view_interactor::compute_holo_views(cgv::render::context& ctx)
{
	double aspect = (double)view_width / view_height;

	if (holo_storage_mode == HM_QUILT) {
		glClearColor(quilt_bg_color.R(), quilt_bg_color.G(), quilt_bg_color.B(), 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		#if COMPUTE == 0
		vi = 0;
		for (quilt_row = 0; quilt_row < quilt_nr_rows; ++quilt_row) {
			for (quilt_col = 0; quilt_col < quilt_nr_cols; ++quilt_col) {

				ivec4 vp(quilt_col * view_width, quilt_row * view_height, view_width, view_height);
				glViewport(vp[0], vp[1], vp[2], vp[3]);
				glScissor(vp[0], vp[1], vp[2], vp[3]);
				glEnable(GL_SCISSOR_TEST);

				current_e = (2.0f * vi) / (nr_holo_views - 1) - 1.0f;
				gl_set_projection_matrix(ctx, current_e, aspect);
				gl_set_modelview_matrix(ctx, current_e, aspect, *this);

				#if EVAL == 0
				// visualise holes
				if (show_holes) 
					mesh->draw_holes(ctx);
				#endif

				switch (multiview_mpx_mode) {
				case MVM_REPROJECT:
					draw_reproject(ctx);
					break;
				case MVM_VWARP_CLOSEST:
					draw_vertex_warp_closest(ctx);
					break;
				case MVM_VWARP:
					draw_vertex_warp(ctx);
					break;
				}
				if (++vi == nr_holo_views)
					break;
			}
			if (vi == nr_holo_views)
				break;
		}
		#endif

		glViewport(0, 0, view_width, view_height);
		glScissor(0, 0, view_width, view_height);
		glDisable(GL_SCISSOR_TEST);
		
		if (multiview_mpx_mode == MVM_COMPUTE || multiview_mpx_mode == MVM_COMPUTE_SPLAT) {
			warp_compute_shader(ctx);
			quilt_resolve_pass_compute_shader(ctx);
		}
	}
	else {
		for (vi = 0; vi < nr_holo_views; ++vi) {
			volume_fbo.attach(ctx, volume_holo_tex, vi, 0, 0);
			glClearColor(0, 0, 0, 1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			#if COMPUTE == 0
			current_e = (2.0f * vi) / (nr_holo_views - 1) - 1.0f;
			gl_set_projection_matrix(ctx, current_e, aspect);
			gl_set_modelview_matrix(ctx, current_e, aspect, *this);

			#if EVAL == 0
			// visualise holes
			if (show_holes)
				mesh->draw_holes(ctx);
			#endif

			switch (multiview_mpx_mode) {
			case MVM_REPROJECT:
				draw_reproject(ctx);
				break;
			case MVM_VWARP_CLOSEST:
				draw_vertex_warp_closest(ctx);
				break;
			case MVM_VWARP:
				draw_vertex_warp(ctx);
				break;
			}
			#endif
		}

		if (multiview_mpx_mode == MVM_COMPUTE || multiview_mpx_mode == MVM_COMPUTE_SPLAT) {
			warp_compute_shader(ctx);
			volume_resolve_pass_compute_shader(ctx);
		}
	}
}

void holo_view_interactor::post_process_surface(cgv::render::context& ctx)
{
	if (generate_hologram) {
		if (!display_fbo.is_created() || display_fbo.get_width() != display_calib.width ||
			display_fbo.get_height() != display_calib.height)
		{
			display_tex.destruct(ctx);
			display_fbo.destruct(ctx);
			display_fbo.create(ctx, display_calib.width, display_calib.height);
			display_tex.create(ctx, cgv::render::TT_2D, display_calib.width, display_calib.height);
			display_fbo.attach(ctx, display_tex);
		}
		cgv::render::shader_program& prog = holo_storage_mode == HM_QUILT ? quilt_prog : volume_prog;
		display_fbo.enable(ctx);
		display_fbo.push_viewport(ctx);
		if (holo_storage_mode == HM_QUILT) {
			prog.set_uniform(ctx, "width", display_calib.width);
			prog.set_uniform(ctx, "nr_views", nr_holo_views);
			prog.set_uniform(ctx, "view_width", view_width);
			prog.set_uniform(ctx, "view_height", view_height);
			prog.set_uniform(ctx, "quilt_nr_cols", quilt_nr_cols);
			prog.set_uniform(ctx, "quilt_nr_rows", quilt_nr_rows);
			prog.set_uniform(ctx, "quilt_width", quilt_width);
			prog.set_uniform(ctx, "quilt_height", quilt_height);
			prog.set_uniform(ctx, "quilt_interpolate", quilt_interpolate);
			quilt_holo_tex.enable(ctx, 0);
			prog.set_uniform(ctx, "quilt_tex", 0);
		}
		else {
			volume_holo_tex.enable(ctx, 0);
			prog.set_uniform(ctx, "volume_tex", 0);
		}
		prog.set_uniform(ctx, "width", display_calib.width);
		prog.set_uniform(ctx, "height", display_calib.height);
		prog.set_uniform(ctx, "length", display_calib.length);
		prog.set_uniform(ctx, "step_x", display_calib.step_x);
		prog.set_uniform(ctx, "step_y", display_calib.step_y);
		prog.set_uniform(ctx, "offset", display_calib.offset);
		prog.set_uniform(ctx, "x_min", display_calib.x_min);
		prog.set_uniform(ctx, "x_max", display_calib.x_max);
		prog.set_uniform(ctx, "y_min", display_calib.y_min);
		prog.set_uniform(ctx, "y_max", display_calib.y_max);
		prog.enable(ctx);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		prog.disable(ctx);
		glDisable(GL_DEPTH_TEST);
		glDepthMask(GL_FALSE);
		glDepthMask(GL_TRUE);
		glEnable(GL_DEPTH_TEST);
		display_fbo.pop_viewport(ctx);
		display_fbo.disable(ctx);
		if (display_write_to_file) {
			display_tex.write_to_file(ctx, "display.png");
			display_write_to_file = false;
			on_set(&display_write_to_file);
		}
		vec4 src_rect(float(blit_offset_x), float(blit_offset_y), float(blit_offset_x + ctx.get_width()),
					  float(blit_offset_y + ctx.get_height()));
		display_fbo.blit_to(ctx, src_rect, cgv::render::BTB_COLOR_BIT, false);
	}

	else {
		if (holo_storage_mode == HM_QUILT) {
			quilt_fbo.blit_to(ctx, BTB_COLOR_BIT, true);
		}
		else {
			volume_fbo.attach(ctx, volume_holo_tex, view_index, 0, 0);
			volume_fbo.blit_to(ctx, BTB_COLOR_BIT, true);
		}
	}

	#if EVAL == 1
	// end time measurement for evaluation
	if (evaluate) {
		glEndQuery(GL_TIME_ELAPSED);

		int done = 0;
		while (!done) {
			glGetQueryObjectiv(time_query, GL_QUERY_RESULT_AVAILABLE, &done);
		}
		glGetQueryObjectui64v(time_query, GL_QUERY_RESULT, &elapsed_time);

		count++;
		accumulated_time += 1000000000.0 / elapsed_time;

		// add the average of 30 frames to the time measurements
		if (count % 30 == 0) {
			time_measurements.push_back(accumulated_time/30.0);
			accumulated_time = 0;
		}

		// rotate by one degree around y-axis
		vec3 new_pos = vec3(eval_pos[0] * cos(count * angle_frac) - eval_pos[2] * sin(count * angle_frac), eval_pos[1],
							eval_pos[2] * cos(count * angle_frac) + eval_pos[0] * sin(count * angle_frac));
		set_eye_keep_extent(new_pos);

		// finish evaluation loop after three full rotations
		if (count >= 3*360) {
			dynamic_cast<fltk_gl_view*>(get_context())->set_void("instant_redraw", "bool", &_off);
			evaluate = false;
			update_member(&evaluate);

			// write all time measurements to .csv file
			std::filesystem::path cwd = std::filesystem::current_path();
			if (cwd.string().find("res") != std::string::npos)
				file.open("../measurements.csv", std::ofstream::in | std::ofstream::app);
			else
				file.open("measurements.csv", std::ofstream::in | std::ofstream::app);
			file << nr_render_views << ", " << multiview_mpx_mode << ", " << mesh->get_number_positions() << ", "
				 << holo_storage_mode;
			for (int i = 0; i < time_measurements.size(); i++) {
				file << ", " << time_measurements[i];
			}
			file << "\n";
			file.close();

			time_measurements.clear();
		}
	}
	#endif
}

void holo_view_interactor::set_up_eval_file() {
	#if EVAL == 1
		int answer = 0;
		if (set_up_file_for_eval) {
			// ask if file should really be flushed to avoid accidental deleting of data	
			answer = cgv::gui::question("Are you sure you want to refresh file?", "Yes, No", 0);
		}

		// if file is not initialised or should be flushed, set up header row
		if (answer == 1 || set_up_file_for_eval == false) {
			std::filesystem::path cwd = std::filesystem::current_path();
			if (cwd.string().find("res") != std::string::npos) 
				file.open("../measurements_new_compute.csv");
			else
				file.open("measurements.csv");
			file << "nr render views, mode, nr vertices, storage";
			for (int i = 0; i < 36; i++) {
				file << ", " << std::to_string(i);
			}
			file << "\n";
			file.close();
		}
	#else
		cgv::gui::message("Please enable EVAL in the code and recompile first.");
	#endif

	set_up_file_for_eval = false;
	update_member(&set_up_file_for_eval);
}

// start or end rotation with instant redraw for evaluation
void holo_view_interactor::toggle_eval()
{
	#if EVAL == 1
		if (evaluate) {
			if (reset_view_for_eval) {
				set_view_dir({0, 0, -1});
				set_view_up_dir({0, 1, 0});
				mesh->focus_mesh();
			}

			eval_pos = get_eye();
			count = 0;
			accumulated_time = 0;
			dynamic_cast<fltk_gl_view*>(get_context())->set_void("instant_redraw", "bool", &_on);
			dynamic_cast<fltk_gl_view*>(get_context())->set_void("vsync", "bool", &_off);
		}

		else {
			dynamic_cast<fltk_gl_view*>(get_context())->set_void("instant_redraw", "bool", &_off);
			evaluate = false;
			update_member(&evaluate);
		}
	#else
		cgv::gui::message("Please enable EVAL in the code and recompile first.");
		evaluate = false;
		update_member(&evaluate);
	#endif
}

///
void holo_view_interactor::draw(cgv::render::context& ctx)
{
	if (show_focus) {
		if (is_viewport_splitting_enabled()) {
			for (unsigned c = 0; c < nr_viewport_columns; ++c) {
				for (unsigned r = 0; r < nr_viewport_rows; ++r) {
					activate_split_viewport(ctx, c, r);
					draw_focus();
					deactivate_split_viewport(ctx);
				}
			}
		}
		else
			draw_focus();
	}
}

void holo_view_interactor::draw_focus()
{
	glLineWidth(1.0f);
	glColor3f(0.5f, 0.5f, 0.5f);
	glBegin(GL_LINES);
	glVertex3dv(get_focus());
	glVertex3dv(get_focus() + dvec3(0.5 * get_y_extent_at_focus(), 0, 0));
	glVertex3dv(get_focus());
	glVertex3dv(get_focus() + dvec3(0, 0.5 * get_y_extent_at_focus(), 0));
	glVertex3dv(get_focus());
	glVertex3dv(get_focus() + dvec3(0, 0, 0.5 * get_y_extent_at_focus()));
	glVertex3dv(get_focus());
	glVertex3dv(get_focus() + dvec3(-0.5 * get_y_extent_at_focus(), 0, 0));
	glVertex3dv(get_focus());
	glVertex3dv(get_focus() + dvec3(0, -0.5 * get_y_extent_at_focus(), 0));
	glVertex3dv(get_focus());
	glVertex3dv(get_focus() + dvec3(0, 0, -0.5 * get_y_extent_at_focus()));
	glEnd();
}



/// return a path in the main menu to select the gui
std::string holo_view_interactor::get_menu_path() const { return "View/Stereo Interactor"; }

/// you must overload this for gui creation
void holo_view_interactor::create_gui()
{
	if (begin_tree_node("View Configuration", zoom_sensitivity, false)) {
		align("\a");
		add_member_control(this, "Use Gamepad", use_gamepad, "toggle");
		add_member_control(this, "Gamepad Emulation", gamepad_emulation, "toggle");
		add_member_control(this, "Pan Sensitivity", pan_sensitivity, "value_slider",
						   "min=0.1;max=10;ticks=true;step=0.01;log=true");
		add_member_control(this, "Zoom Sensitivity", zoom_sensitivity, "value_slider",
						   "min=0.1;max=10;ticks=true;step=0.01;log=true");
		add_member_control(this, "Rotate Sensitivity", rotate_sensitivity, "value_slider",
						   "min=0.1;max=10;ticks=true;step=0.01;log=true");
		add_member_control(this, "Deadzone", deadzone, "value_slider", "min=0;max=1;ticks=true;step=0.0001;log=true");
		add_member_control(this, "Viewport Shrinkage", viewport_shrinkage, "value_slider", "min=0;max=40;ticks=true");
		add_member_control(this, "Show Focus", show_focus, "check");
		align("\b");
		end_tree_node(zoom_sensitivity);
	}
	if (begin_tree_node("Holo Parameters", multiview_mpx_mode, true)) {
		align("\a");
		if (begin_tree_node("Stereo", eye_distance, false)) {
			align("\a");
			add_member_control(this, "Eye Distance", eye_distance, "value_slider",
							   "min=0.01;max=5;ticks=true;step=0.00001;log=true");
			add_member_control(this, "Parallax Zero Scale", parallax_zero_scale, "value_slider",
							   "min=0.03;max=1;ticks=true;step=0.001;log=true");
			add_member_control(this, "Stereo Translate in Model View", stereo_translate_in_model_view, "check");
			align("\b");
			end_tree_node(eye_distance);
		}
		if (begin_tree_node("Display", display_calib, false)) {
			align("\a");
			add_member_control(this, "Width", display_calib.width, "value_slider", "min=1920;max=8192;ticks=true");
			add_member_control(this, "Height", display_calib.height, "value_slider", "min=1080;max=4096;ticks=true");
			add_gui("X_min", display_calib.x_min, "", "options='min=0;max=0.1;ticks=true'");
			add_gui("X_max", display_calib.x_max, "", "options='min=0.9;max=1;ticks=true'");
			add_gui("Y_min", display_calib.y_min, "", "options='min=0;max=0.1;ticks=true'");
			add_gui("Y_max", display_calib.y_max, "", "options='min=0.9;max=1;ticks=true'");
			add_member_control(this, "Length", display_calib.length, "value_slider",
							   "min=42;max=43;ticks=true;step=0.000001");
			add_member_control(this, "Step_x", display_calib.step_x, "value_slider",
							   "min=35;max=36;ticks=true;step=0.000001");
			add_member_control(this, "Step_y", display_calib.step_y, "value_slider",
							   "min=0.9;max=1.1;ticks=true;step=0.000001");
			add_gui("Offset", display_calib.offset, "", "options='min=37;max=43;ticks=true'");
			end_tree_node(display_calib);
		}
		if (begin_tree_node("Rendering", multiview_mpx_mode, true)) {
			align("\a");
			add_member_control(this, "Use Optimal", choose_optimal, "toggle");
			add_member_control(this, "Render Multiplexing", multiview_mpx_mode, "dropdown",
							   "enums='single view, conventional, geometry, reproject, vertex warp, vertex warp closest, "
							   "compute warp splat, compute warp'");
			add_member_control(this, "Storage Mode", holo_storage_mode, "dropdown",
							   "enums='single view,quilt,volume'");
			add_member_control(this, "View Width", view_width, "value_slider", "min=640;max=2000;ticks=true");
			add_member_control(this, "View Height", view_height, "value_slider", "min=480;max=1000;ticks=true");
			add_member_control(this, "Number Rendered Views", nr_render_views, "value_slider",
							   "min=1;max=3;ticks=true");
			add_member_control(this, "Number Hologram Views", nr_holo_views, "value_slider",
							   "min=2;max=100;ticks=true");
			add_member_control(this, "Epsilon", epsilon, "value_slider", "min=0;max=0.1;step=0.00001;ticks=true");
			add_member_control(this, "View Index", view_index, "value_slider", "min=0;max=44;ticks=true");
			add_member_control(this, "Blit Offset x", blit_offset_x, "value_slider", "min=0;max=1000;ticks=true");
			add_member_control(this, "Blit Offset y", blit_offset_y, "value_slider", "min=0;max=1000;ticks=true");
			add_member_control(this, "Discard Rubber-Sheets", dis_artefacts, "check");
			add_member_control(this, "Show Holes", show_holes, "check");
			add_member_control(this, "Reset View When Evaluating", reset_view_for_eval, "check");
			connect_copy(add_member_control(this, "Set Up File", set_up_file_for_eval, "toggle")->value_change,
						 rebind(this, &holo_view_interactor::set_up_eval_file));
			connect_copy(add_member_control(this, "Start Evaluation Run", evaluate, "toggle")->value_change,
						 rebind(this, &holo_view_interactor::toggle_eval));
			add_member_control(this, "Generate Hologram", generate_hologram, "toggle");
			add_member_control(this, "Write To File", display_write_to_file, "toggle");
			end_tree_node(multiview_mpx_mode);
		}
		if (begin_tree_node("Quilt", quilt_bg_color, true)) {
			align("\a");
			add_member_control(this, "Background", quilt_bg_color);
			add_member_control(this, "Use Offline Texture", quilt_use_offline_texture, "toggle");
			add_member_control(this, "Columns", quilt_nr_cols, "value_slider", "min=0;max=10;ticks=true");
			add_member_control(this, "Rows", quilt_nr_rows, "value_slider", "min=0;max=20;ticks=true");
			add_member_control(this, "Interpolate", quilt_interpolate, "check");
			add_member_control(this, "Write To File", quilt_write_to_file, "toggle");
			align("\b");
			end_tree_node(quilt_bg_color);
		}
		if (begin_tree_node("Volume", true, true)) {
			align("\a");
			add_member_control(this, "Write To File", volume_write_to_file, "toggle");
			align("\b");
			end_tree_node(quilt_bg_color);
		}
		align("\b");
		end_tree_node(multiview_mpx_mode);
	}
	if (begin_tree_node("Current View", view::focus(0), true)) {
		align("\a");
		add_member_control(this, "Focus x", view::focus(0), "value_input", "w=50;min=-10;max=10;ticks=true", " ");
		add_member_control(this, "y", view::focus(1), "value_input", "w=50;min=-10;max=10;ticks=true", " ");
		add_member_control(this, "z", view::focus(2), "value_input", "w=50;min=-10;max=10;ticks=true");
		///
		add_gui("View Dir", view_dir, "direction", "options='min=-1;max=1;ticks=true'");
		add_gui("View Up Dir", view_up_dir, "direction", "options='min=-1;max=1;ticks=true'");
		connect_copy(add_control("Fix View Up Dir", fix_view_up_dir, "check")->value_change,
					 rebind(this, &holo_view_interactor::on_rotation_change));

		add_member_control(this, "y View Angle", y_view_angle, "value_slider", "min=0;max=90;ticks=true;log=true");
		add_member_control(this, "y Extent at Focus", y_extent_at_focus, "value_slider",
						   "min=0;max=100;ticks=true;log=true;step=0.0001");
		add_member_control(this, "z Near", z_near, "value_slider", "min=0;max=100;log=true;step=0.00001");
		add_member_control(this, "z Far", z_far, "value_slider", "min=0;max=10000;log=true;step=0.00001");
		add_member_control(this, "Clip Relative To Extent", clip_relative_to_extent, "check");
		align("\b");
		end_tree_node(view::focus(0));
	}
}

void holo_view_interactor::on_set(void* m)
{
	if (m == &nr_holo_views) {
		if (view_index >= nr_holo_views)
			view_index = nr_holo_views - 1;
		if (find_control(view_index))
			find_control(view_index)->set("max", nr_holo_views - 1);
		update_member(&view_index);
		change_nr_holo_views = true;
	}
	// make sure MVM_SINGLE is not shown as quilt or volume
	if (multiview_mpx_mode == MVM_SINGLE && holo_storage_mode != HM_SINGLE)
		holo_storage_mode = HM_SINGLE;
	if (m == &multiview_mpx_mode && (multiview_mpx_mode == MVM_COMPUTE || multiview_mpx_mode == MVM_COMPUTE_SPLAT))
		update_defines = true;
	// change size of shader storage buffer in case viewport changes
	if (m == &view_width || m == &view_height || m == &nr_holo_views || m == &quilt_nr_cols) {
		glDeleteBuffers(1, &ssbo);
		glGenBuffers(1, &ssbo);
		glNamedBufferData(ssbo,
						  GLsizeiptr(sizeof(unsigned int) * view_width * view_height *
									 std::ceil(nr_holo_views / (float)quilt_nr_cols) * quilt_nr_cols),
							nullptr, GL_DYNAMIC_COPY);
	}
	// in case the optimal solution is chosen, the conventional approach should be used for smaller scenes. The compute approach with splat and three source views
	// should be chosen for larger scenes. The volume mode is used for the storage modus because it is slightly faster for the conventional approach for smaller scenes and
	// slightly faster for the compute splat approach for larger scenes
	if (choose_optimal) {
		int num_vertices = mesh->get_number_positions();
		if (num_vertices < 750000) {
			multiview_mpx_mode = MVM_CONVENTIONAL;
		}
		else {
			multiview_mpx_mode = MVM_COMPUTE_SPLAT;
			nr_render_views = 3;
			update_member(&nr_render_views);
		}
		update_member(&multiview_mpx_mode);
		holo_storage_mode = HM_VOLUME;
		update_member(&holo_storage_mode);
	}
	update_member(m);
	post_redraw();
}

void holo_view_interactor::set_z_near(double z)
{
	cgv::render::clipped_view::set_z_near(z);
	update_member(&z_near);
	post_redraw();
}
void holo_view_interactor::set_z_far(double z)
{
	cgv::render::clipped_view::set_z_far(z);
	update_member(&z_far);
	post_redraw();
}
void holo_view_interactor::set_default_view()
{
	cgv::render::clipped_view::set_default_view();
	for (unsigned c = 0; c < 3; ++c) {
		update_member(&view_dir[c]);
		update_member(&view_up_dir[c]);
		update_member(&cgv::render::view::focus[c]);
	}
	update_member(&y_extent_at_focus);
	post_redraw();
}

/// you must overload this for gui creation
bool holo_view_interactor::self_reflect(cgv::reflect::reflection_handler& srh)
{
	return srh.reflect_member("enable_messages", enable_messages) && srh.reflect_member("use_gamepad", use_gamepad) &&
		   srh.reflect_member("gamepad_emulation", gamepad_emulation) && srh.reflect_member("deadzone", deadzone) &&
		   srh.reflect_member("pan_sensitivity", pan_sensitivity) &&
		   srh.reflect_member("rotate_sensitivity", rotate_sensitivity) &&
		   srh.reflect_member("zoom_sensitivity", zoom_sensitivity) && srh.reflect_member("focus_x", view::focus(0)) &&
		   srh.reflect_member("focus_y", view::focus(1)) && srh.reflect_member("focus_z", view::focus(2)) &&
		   srh.reflect_member("focus", view::focus) && srh.reflect_member("eye_distance", eye_distance) &&
		   srh.reflect_member("parallax_zero_scale", parallax_zero_scale) &&
		   srh.reflect_member("stereo_translate_in_model_view", stereo_translate_in_model_view) &&
		   srh.reflect_member("view_dir_x", view_dir(0)) && srh.reflect_member("view_dir_y", view_dir(1)) &&
		   srh.reflect_member("view_dir_z", view_dir(2)) && srh.reflect_member("view_dir", view_dir) &&
		   srh.reflect_member("fix_view_up_dir", fix_view_up_dir) && srh.reflect_member("up_dir_x", view_up_dir(0)) &&
		   srh.reflect_member("up_dir_y", view_up_dir(1)) && srh.reflect_member("up_dir_z", view_up_dir(2)) &&
		   srh.reflect_member("up_dir", view_up_dir) && srh.reflect_member("y_view_angle", y_view_angle) &&
		   srh.reflect_member("extent", y_extent_at_focus) && srh.reflect_member("z_near", z_near) &&
		   srh.reflect_member("z_far", z_far) && srh.reflect_member("two_d_enabled", two_d_enabled) &&
		   srh.reflect_member("show_focus", show_focus) &&
		   srh.reflect_member("clip_relative_to_extent", clip_relative_to_extent) &&

		   srh.reflect_member("display_width", display_calib.width) &&
		   srh.reflect_member("display_height", display_calib.height) &&
		   srh.reflect_member("display_x_min", display_calib.x_min) &&
		   srh.reflect_member("display_x_max", display_calib.x_max) &&
		   srh.reflect_member("display_y_min", display_calib.y_min) &&
		   srh.reflect_member("display_y_max", display_calib.y_max) &&
		   srh.reflect_member("display_length", display_calib.length) &&
		   srh.reflect_member("display_step_x", display_calib.step_x) &&
		   srh.reflect_member("display_step_y", display_calib.step_y) &&
		   srh.reflect_member("display_offset", display_calib.offset) &&
		   srh.reflect_member("multiview_mpx_mode", multiview_mpx_mode) &&
		   srh.reflect_member("holo_mpx_mode", holo_storage_mode) && srh.reflect_member("view_width", view_width) &&
		   srh.reflect_member("view_height", view_height) && srh.reflect_member("nr_render_views", nr_render_views) &&
		   srh.reflect_member("nr_holo_views", nr_holo_views) && srh.reflect_member("view_index", view_index) &&
		   srh.reflect_member("blit_offset_x", blit_offset_x) && srh.reflect_member("blit_offset_y", blit_offset_y) &&
		   srh.reflect_member("generate_hologram", generate_hologram) &&
		   srh.reflect_member("quilt_bg_color", quilt_bg_color) &&
		   srh.reflect_member("quilt_use_offline_texture", quilt_use_offline_texture) &&
		   srh.reflect_member("quilt_nr_cols", quilt_nr_cols) && srh.reflect_member("quilt_nr_rows", quilt_nr_rows) &&
		   srh.reflect_member("quilt_interpolate", quilt_interpolate) &&
		   srh.reflect_member("dis_artefacts", dis_artefacts) && srh.reflect_member("show_holes", show_holes) &&
		   srh.reflect_member("epsilon", epsilon);
}

// register custom holo-view interactor
#include <cgv/base/register.h>
cgv::base::object_registration_1<holo_view_interactor, const char*> obj1("Holo Interactor",
																		 "registration of holo interactor");
