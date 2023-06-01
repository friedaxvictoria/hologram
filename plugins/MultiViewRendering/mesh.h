#pragma once

#include <cgv/media/mesh/simple_mesh.h>
#include <cgv/base/base.h>
#include <cgv/render/drawable.h>
#include <cgv/gui/provider.h>
#include <cgv_gl/sphere_renderer.h>
#include <cgv_gl/cone_renderer.h>
#include <cgv_gl/gl/mesh_render_info.h>
#include <cgv_gl/box_wire_render_data.h>
#include <cgv/render/managed_frame_buffer.h>
#include <cgv/render/stereo_view.h>
#include <cgv/gui/event_handler.h>

#include "lib_begin.h"
#include "tessellator.h"

using namespace cgv::base;
using namespace cgv::signal;
using namespace cgv::gui;
using namespace cgv::math;
using namespace cgv::render;
using namespace cgv::utils;
using namespace cgv::media::illum;

class mesh_viewer : 
									 public node,
									 public drawable,
									 public provider,
									 public event_handler
   {
  public:
	typedef cgv::math::fvec<float, 3> vec3;
	typedef cgv::math::fvec<double, 3> dvec3;
	typedef cgv::math::fmat<double, 3, 3> dmat3;
	typedef cgv::math::fmat<double, 4, 4> dmat4;
	typedef cgv::media::mesh::simple_mesh<float> mesh_type;
	typedef mesh_type::idx_type idx_type;
	typedef mesh_type::vec3i vec3i;

  protected:
	float fac = 0.00001;
	std::string mesh_filename;
	mesh_type M;
	cgv::render::mesh_render_info mesh_info, mesh_for_geo_info, mesh_for_holes_info;
	cgv::render::box3 M_bbox;
	cgv::render::box_wire_render_data<> M_bbox_rd;
	bool meshfile_supplies_colors, invent_missing_colors = false;

	bool update_view_after_mesh_processed = false;

	bool show_bbox = false;

	bool show_surface = true;
	CullingMode cull_mode;
	ColorMapping color_mapping;
	rgb surface_color;
	IlluminationMode illumination_mode;

	bool show_vertices = false;
	sphere_render_style sphere_style;
	attribute_array_manager sphere_aam;

	bool show_wireframe = false;
	cone_render_style cone_style;
	attribute_array_manager cone_aam;

	cgv::render::stereo_view* view = nullptr;

	vec3 cam_dir;

	shader_program geometry_shader; // the mesh-based baseline approach
	shader_program holes_shader;	// for displaying the holes

	// variables for geometry shader render call
	float zero_parallax, eye_distance, num_holo_views, eye;
	bool with_geometry = false;

	public:
		mesh_viewer();
		void on_register();
		std::string get_type_name() const;
		void process_mesh_for_rendering(context& ctx, bool update_view);
		void ensure_mesh_colors();
		bool handle(cgv::gui::event& e);
		void stream_help(std::ostream& os);
		void on_set(void* member_ptr);
		void create_conway_polyhedron();
		bool init(context& ctx);
		void clear(context& ctx);
		void init_frame(context& ctx);
		void draw_holes(context& ctx);
		void draw_geometry_shader(context& ctx);
		void set_params_for_gemoetry(float zero_parallax, float eye_distance, float eye, float num_holo_views);
		void draw_surface(context& ctx, bool opaque_part);
		void draw(context& ctx);
		void finish_frame(context& ctx);
		bool self_reflect(cgv::reflect::reflection_handler& srh);
		void create_gui();
		void focus_mesh();
   };
