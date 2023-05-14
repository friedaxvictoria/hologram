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

	float eye_distance;

	vec3 cam_dir;

	////
	// 3D Image Warp baseline testing fields

	struct
	{
		enum WarpMode { BASELINE, IMAGE_WARP };
		WarpMode warp_mode = IMAGE_WARP;

		// Render targets for each fully rendered view - we allocate enough for 3 viewpoints:
		// 0-left, 1-center, 2-right
		cgv::render::managed_frame_buffer render_fbo[3];

		// our base projection matrix
		mat4 inv_mat_proj_render[3];
		mat4 proj_for_render;

		// image warping test shaders
		shader_program baseline_shader; // the mesh-based baseline approach
		shader_program holes_shader;	// for displaying the holes
		shader_program warp_shader;		// the image warping approach
		shader_program geometry_shader;		// the image warping approach

		// image warping shader parameters
		bool prune_heightmap =
			  true; // discard heightmap fragments which don't represent valid geometry (from "empty" areas)
		float heightmap_oversampling =
			  2.f; // oversampling factor of the heightmap (to be able to resolve finer details)

		// Mesh for the heightmap geometry (actually the same for all three views,
		// since the topology never changes!)
		GPUgeometry heightmap_baseline, heightmap_warp;

		// transformation matrix for positioning the heightmap at the plane behind the scene from
		// which it was shot
		mat4 heightmap_trans;
		mat4 modelview_source;

		// whether to render the heightmap
		bool render_heightmap = false;

		// indicates that a snapshot of the current view should be safed in the heightmap
		bool shoot_heightmap;

		int visible_view = 22;

		int nr_renders = 3;
		float render_offset[3];

		bool show_holes = true;
		bool with_interpolated_holes = false;
		bool nr_rendered_views_changed = false;

		float epsilon = 0.02;

		bool ortho = false;

		float x_ext, y_ext, znear;
		vec3 eye_source[3], eye_target;

		mat3 p_1[3];
	} test;

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
		void draw_geometry_shader(context& ctx, float zero_parallax, float eye_distance, int eye,
											   float eye_offset);
		void draw_surface(context& ctx, bool opaque_part);
		void draw(context& ctx);
		void finish_frame(context& ctx);
		bool self_reflect(cgv::reflect::reflection_handler& srh);
		void on_shoot(void);
		void nr_rendered_views_change(void);
		void create_gui();
   };
