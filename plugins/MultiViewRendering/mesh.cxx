#include <cgv/base/node.h>
#include <cgv/defines/quote.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/drawable.h>
#include <cgv/render/clipped_view.h>
#include <cgv_gl/gl/gl.h>
#include <cgv_gl/sphere_renderer.h>
#include <cgv_gl/cone_renderer.h>
#include <cgv_gl/gl/mesh_render_info.h>
#include <cgv_gl/gl/gl_context.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/application.h>
#include <cgv/media/color_scale.h>
#include <cgv/media/mesh/simple_mesh.h>
#include <cgv/gui/key_event.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/dialog.h>
#include <cgv/math/ftransform.h>
#include <cgv/render/render_buffer.h>
#include <cgv/render/frame_buffer.h>


using namespace cgv::base;
using namespace cgv::signal;
using namespace cgv::gui;
using namespace cgv::math;
using namespace cgv::render;
using namespace cgv::utils;
using namespace cgv::media::illum;

/*
This example illustrates how to use cgv::mesh::simple_mesh<T> and cgv::render::mesh_render_info
to render a polygonal mesh in potentially textured surface rendering mode and in wireframe mode.
The mesh_render_info manages vertex buffers for position, normal, tex_coord and color attributes
and an element index buffer for triangle based surface and edge based wireframe rendering. 
Furthermore, an attribute_array_binding object is managed.

The typical workflow is

PREPARATION (implemented in init or init_frame method of drawable)
- construct simple_mesh data structure with new_position(), new_normal(), new_tex_coord(), 
  start_face() and new_corner() methods. If also per vertex colors are needed, allocate
  them with dynamically choosable type through the ensure_colors() method of the base
  class cgv::media::colored_model and the same number of colors as you have positions.
  Then individual colors can be set with the set_color() functions and retreived later on
  with the put_color() functions. 
  An alternative is to read a simple mesh from an obj-file with the read() method.
  Optionally, you can compute vertex normals with the compute_vertex_normals() method after
  construction of positions, faces and corners.
  [this step can be implemented also outside of the drawable methods as no render context is needed]
- construct the vertex buffers of the mesh_render_info object from the simple_mesh with its 
  construct_vbos() method
- choose a shader program and bind the attribute_array_binding of the mesh_render_info struct to
  the attribute locations of the program with the bind() method. This will also bind the element
  index buffer.

RENDERING PHASE (implemented in draw and finish_draw method)
- configure uniforms and constant vertex attributes of shader program 
- call render_mesh() method of mesh_render_info and pass chosen shader program to it.
  This will enable and disable the shader program as needed. render_mesh() has two optional
  boolean parameters that allow to ignore transparent or opaque mesh parts during rendering.
  If the mesh has transparent parts, one typically renders first all opaque parts and in a 
  second pass or in the finish_draw method all transparent parts. This will not produce a 
  correct visibility ordering of the transparent parts but will blend over the transparent
  parts over the opaque parts. If not done like this, transparent parts can occlude opaque
  parts due to the z-buffer algorithm used for visibility sorting.

The mesh_render_info class provides part based rendering methods for meshes with several
materials and several groups as supported by the obj file format with the render_mesh_part()
method. Furthermore, there are two draw methods draw_surface() and draw_wireframe(). Both 
methods do not enable any shader program and assume the vertex attribute locations of the
used shader program to be the same as of the shader program passed to the bind() method in
the preparation stage. In this example the draw_wireframe() method is used for wireframe 
rendering.

The example furthermore supports reading of obj files and can serve as simple mesh viewer.
Some image formats used for 3d models are not support by the cmi_io plugin used by default
in the examples plugin. More image formats are provided by the the cmi_devIL plugin which is
only part of the cgv_support project tree available on demand.
*/

enum WarpingMode { SIMPLE_INTERPOLATION, SIMPLE_WARPING };

class mesh_viewer : public node, public drawable, public provider
{
  public:
	typedef cgv::media::mesh::simple_mesh<float> mesh_type;
	typedef mesh_type::idx_type idx_type;
	typedef mesh_type::vec3i vec3i;

  protected:
	clipped_view* view_ptr;

	std::string mesh_filename;
	mesh_type M;
	cgv::render::mesh_render_info mesh_info;
	bool meshfile_supplies_colors, invent_missing_colors = false;

	bool update_view_after_mesh_processed = false;

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

	////
	// 3D Image Warp baseline testing fields

	// Render targets for each fully rendered view
	cgv::render::texture render_color[3];        // we allocate enough for
	cgv::render::render_buffer render_depth[3];  // 3 viewpoints: 0-left,
	cgv::render::frame_buffer render_fbo[3];     // 1-center, 2-right

	// Mesh for the heightmap geometry (actually the same for all three views,
	// since the topology never changes!)
	mesh_type heightmap;
	cgv::render::mesh_render_info heightmap_info;

	cgv::render::shader_program baseline_prog;

	vec3 render_eye_position;
	double y_extend;
	int visible_view = 22;
	WarpingMode warping_mode;	

public:
	/// the constructor
	mesh_viewer() : node("mesh_viewer"), render_depth{"[D]", "[D]", "[D]"}
	{
		cull_mode = CM_BACKFACE;
		color_mapping = cgv::render::CM_COLOR;
		surface_color = rgb(0.75f, 0.25f, 1.0f);
		illumination_mode = IM_ONE_SIDED;
		warping_mode = SIMPLE_INTERPOLATION;

		sphere_style.surface_color = rgb(0.8f, 0.3f, 0.3f);
		cone_style.surface_color = rgb(0.6f, 0.5f, 0.4f);
	}

	/// reflect the name of our class
	std::string get_type_name() const {
		return "mesh_viewer";
	}

	/// called when an instance of this class is registered with the Framework
	void on_register()
	{
		// --NOTE-- obviously, a window can only have one title, so it should be set in the most central, "main"
		//          class that the project has (which this one, as it just provides some geometry to display, might
		//          not be anymore at some point)
		cgv::gui::application::get_window(0)->set("title", "Multi-View Rendering for Holographic Displays");
	}

	/// helper function that acts as a single point where processing of the mesh into a renderable form is happening
	void process_mesh_for_rendering(context& ctx, bool update_view = true)
	{
		// adjust size of vertex and edge glyphs to loaded mesh
		sphere_style.radius = float(0.05 * sqrt(M.compute_box().get_extent().sqr_length() / M.get_nr_positions()));
		on_set(&sphere_style.radius);
		cone_style.radius = 0.5f * sphere_style.radius;
		on_set(&cone_style.radius);

		// if requested by the user, invent per-vertex colors in case the mesh has none
		meshfile_supplies_colors = M.has_colors();
		if (invent_missing_colors)
			ensure_mesh_colors();
		// compute mesh normals if not yet present
		if (!M.has_normals())
			M.compute_vertex_normals();

		// [re-]compute mesh render info
		mesh_info.destruct(ctx);
		mesh_info.construct(ctx, M);
		// bind mesh attributes to standard surface shader program
		mesh_info.bind(ctx, ctx.ref_surface_shader_program(true), true);
		mesh_info.bind_wireframe(ctx, ref_cone_renderer(ctx).ref_prog(), true);

		// update sphere attribute array manager
		sphere_renderer& sr = ref_sphere_renderer(ctx);
		sr.enable_attribute_array_manager(ctx, sphere_aam);
		sr.set_position_array(ctx, M.get_positions());
		if (M.has_colors())
			sr.set_color_array(ctx, *reinterpret_cast<const std::vector<rgb>*>(M.get_color_data_vector_ptr()));

		// adjust camera parameters when requested
		update_view_after_mesh_processed |= update_view;

		// tessellation of heightmap
		// add positions to mesh - one position for each pixel
		for (int x = 0; x < ctx.get_width(); ++x) {
			for (int y = 0; y < ctx.get_height(); ++y) {
				heightmap.new_position(vec3((float const&)x, (float const&)y, 0));
			}
		}

		// create indices for a quad mesh
		std::vector<int> indices;
		for (int x = 0; x < ctx.get_width() - 1; ++x) {
			for (int y = 0; y < ctx.get_height()-1; ++y) {
				indices.push_back(x * (ctx.get_height()+1) + y);
				indices.push_back((x + 1) * (ctx.get_height()+1) + y);
				indices.push_back(x * (ctx.get_height() + 1) + y+1);
				indices.push_back((x + 1) * (ctx.get_height() + 1) + y+1);
			}
		}

		// create faces for the quad mesh based off of the indices
		int i = 0;
		for (int fi = 0; fi < (ctx.get_width() - 1) * (ctx.get_height()-1); ++fi) {
			heightmap.start_face();
			for (int ci = 0; ci < 4; ++ci) {
				heightmap.new_corner(indices[i], fi);
				i++;
			}
		}

		// compute heightmap and bind it to the baseline program
		/* heightmap_info.destruct(ctx);
		heightmap_info.construct(ctx, heightmap);
		heightmap_info.bind(ctx, baseline_prog, true);*/
		

	


		// ensure that materials are presented in gui
		post_recreate_gui();
	}

	/// helper function that will make sure we have a per-vertex color attribute
	void ensure_mesh_colors()
	{
		if (M.has_colors())
			return;
		const int nr_positions = M.get_nr_positions();
		M.ensure_colors(cgv::media::CT_RGB, nr_positions);
		double dummy;
		#pragma omp for
		for (int i=0; i<nr_positions; i++) {
			M.set_color(i, cgv::media::color_scale(modf(double(20*i)/double(nr_positions-1), &dummy)));
		}
	}

	/// react to our class fields being set via the GUI or via reflection (e.g. from a config file)
	void on_set(void *member_ptr)
	{
		if (member_ptr == &mesh_filename)
		{
			mesh_type tmp;
			if (tmp.read(mesh_filename)) {
				M = std::move(tmp);
				process_mesh_for_rendering(*get_context());
			}
		}
		if (member_ptr == &invent_missing_colors) {
			if (!meshfile_supplies_colors && !invent_missing_colors)
				M.destruct_colors();
			process_mesh_for_rendering(*get_context(), false);
		}

		post_redraw();
		update_member(member_ptr);
	}

	/// clears whatever mesh is currently loaded an creates a Conway polyhedron instead
	void create_conway_polyhedron ()
	{
		M.clear();
		M.construct_conway_polyhedron("dtI");
		mesh_filename = "<generated>";
		process_mesh_for_rendering(*get_context());
	}


	/// perform initialization that can not be done in the constructor as it requires a fully functional graphics
	/// context
	bool init(context &ctx)
	{
		// init render components
		ref_sphere_renderer(ctx, 1);
		ref_cone_renderer(ctx, 1);
		sphere_aam.init(ctx);
		cone_aam.init(ctx);

		// in the blank state (without anything loaded), we just display a simple Conway polyhedron
		create_conway_polyhedron();

		if (!baseline_prog.build_program(ctx, "baseline.glpr", true))
			return false;

		// repost success
		return true;
	}

	/// unload any currently loaded mesh data
	void clear(context &ctx)
	{
		ref_cone_renderer(ctx, -1);
		ref_sphere_renderer(ctx, -1);
		sphere_aam.destruct(ctx);
		cone_aam.destruct(ctx);
	}

	/// perform any kind of operation that should take place before all drawables start executing their ::draw()
	/// methods
	void init_frame(context &ctx)
	{
		////
		// SECTION: 3D image warping baseline test

		// convenience shortcuts
		const unsigned width = ctx.get_width(), height = ctx.get_height();

		// Check if render fbo needs re-initialization (we always change all 3 at the same time, so we just check the center one)
		if (!render_fbo[1].is_created() || render_fbo[1].get_width() != width || render_fbo[1].get_height() != height)
		{
			// (re-)initialize all three fbos and associated resources
			for (unsigned i=0; i<3; i++)
			{
				render_fbo[i].destruct(ctx);
				render_color[i].destruct(ctx);
				render_depth[i].destruct(ctx);
				render_color[i].create(ctx, TT_2D, width, height);
				render_depth[i].create(ctx, width, height);
				render_fbo[i].create(ctx, width, height);
				render_fbo[i].attach(ctx, render_color[i]);
				render_fbo[i].attach(ctx, render_depth[i]);
			}
		}

		// END: 3D image warping baseline test
		////


		if (update_view_after_mesh_processed)
		{
			// focus view on new mesh
			view_ptr = dynamic_cast<clipped_view*>(find_view_as_node());

			if (view_ptr)
			{
				box3 box = M.compute_box();
				view_ptr->set_scene_extent(box);
				view_ptr->set_focus(box.get_center());
				view_ptr->set_y_extent_at_focus(box.get_extent().length());
			}
			update_view_after_mesh_processed = false;
		}
	}

	/// draw the mesh surface
	void draw_surface(context &ctx, bool opaque_part)
	{
		// remember current culling setting
		GLboolean is_culling = glIsEnabled(GL_CULL_FACE);
		GLint cull_face;
		glGetIntegerv(GL_CULL_FACE_MODE, &cull_face);

		// ensure that opengl culling is identical to shader program based culling
		if (cull_mode > 0) {
			glEnable(GL_CULL_FACE);
			glCullFace(cull_mode == CM_BACKFACE ? GL_BACK : GL_FRONT);
		}
		else
			glDisable(GL_CULL_FACE);

		// choose a shader program and configure it based on current settings
		shader_program &prog = ctx.ref_surface_shader_program(true);
		prog.set_uniform(ctx, "culling_mode", (int)cull_mode);
		prog.set_uniform(ctx, "map_color_to_material", (int)color_mapping);
		prog.set_uniform(ctx, "illumination_mode", (int)illumination_mode);
		// set default surface color for color mapping which only affects 
		// rendering if mesh does not have per vertex colors and color_mapping is on
		//prog.set_attribute(ctx, prog.get_color_index(), surface_color);
		ctx.set_color(surface_color);
		// render the mesh from the vertex buffers with selected program
		mesh_info.draw_all(ctx, opaque_part, !opaque_part);

		// recover opengl culling mode
		if (is_culling)
			glEnable(GL_CULL_FACE);
		else
			glDisable(GL_CULL_FACE);
		glCullFace(cull_face);
	}

	/// the "main" draw method
	void draw(context &ctx)
	{
		////
		// SECTION: 3D image warping baseline test

		// --NOTE-- bind appropriate render_fbo here (e.g. just render_fbo[1] for the
		//          initial mono depth map test case)

		// END: 3D image warping baseline test
		////

		render_fbo[1].enable(ctx, 0);
		render_fbo[1].push_viewport(ctx);

		render_eye_position = view_ptr->get_eye();
		y_extend = view_ptr->get_y_extent_at_depth(render_eye_position[0], true);

		float aspect = (visible_view / 22.0 - 1.0);
		float offset_for_current_view = y_extend * aspect;
		vec3 current_eye_position =
			  vec3(render_eye_position[0] + offset_for_current_view, render_eye_position[1], render_eye_position[2]);

		auto MVPW_source = ctx.get_modelview_projection_window_matrix(), 
			MVP_source = ctx.get_modelview_matrix() * ctx.get_projection_matrix();
		view_ptr->set_eye_keep_extent(vec3(render_eye_position[0] + offset_for_current_view, render_eye_position[1],
									  render_eye_position[2]));
		auto MVPW_target = ctx.get_modelview_projection_window_matrix(),
			 MVP_target = ctx.get_modelview_matrix() * ctx.get_projection_matrix();
		view_ptr->set_eye_keep_extent(render_eye_position);

		if (show_vertices)
		{
			sphere_renderer &sr = ref_sphere_renderer(ctx);
			sr.set_render_style(sphere_style);
			sr.enable_attribute_array_manager(ctx, sphere_aam);
			sr.render(ctx, 0, M.get_nr_positions());
			sr.disable_attribute_array_manager(ctx, sphere_aam);
		}
		if (show_wireframe)
		{
			cone_renderer &cr = ref_cone_renderer(ctx);
			cr.set_render_style(cone_style);
			if (cr.enable(ctx)) {
				mesh_info.draw_wireframe(ctx);
				cr.disable(ctx);
			}
		}
		if (show_surface)
			draw_surface(ctx, true);

		
		////
		// SECTION: 3D image warping baseline test

		// --NOTE-- Disable the render_fbo to make the "main" framebuffer active again.
		//          Render a single test view (e.g. selectable via GUI) by rendering the
		//          heightmap mesh using the corresponding transformation while the color
		//          and depth buffer for the central view are actively bound to the appropriate
		//          texture units (so the vertex shader can offset and color the vertices
		//          accordingly).
		//          For the other warping approaches, you would only draw a screen quad (see
		//          holo_view_interactor!) and the fragment shader would do the main work

		// END: 3D image warping baseline test
		////

		render_fbo[1].disable(ctx);

		/* shader_program& prog = baseline_prog;
		prog.set_attribute(ctx, 0, heightmap.get_positions());
		prog.set_uniform(ctx, "mvpw_source", MVPW_source);
		prog.set_uniform(ctx, "mvpw_target", MVPW_target);
		prog.set_uniform(ctx, "mvp_source", MVP_source);
		prog.set_uniform(ctx, "mvp_target", MVP_target);
		prog.set_uniform(ctx, "width", (float)ctx.get_width());
		prog.set_uniform(ctx, "height", (float)ctx.get_height());
		prog.set_uniform(ctx, "eye_pos_rendered", render_eye_position);
		prog.set_uniform(ctx, "eye_pos_current", current_eye_position);
		prog.set_uniform(ctx, "current_view", visible_view);
		prog.set_uniform(ctx, "depth_tex", render_depth[1]);
		prog.set_uniform(ctx, "colour_tex", render_color[1]);
		prog.set_uniform(ctx, "warping_mode", (int)warping_mode);

		prog.enable(ctx);
		heightmap_info.draw_all(ctx);
		prog.disable(ctx);
		*/
	}

	/// perform any kind of operation that should take place after all drawables have executed their ::draw()
	/// methods
	void finish_frame(context &ctx)
	{
		if (show_surface)
			draw_surface(ctx, false);
	}

	/// reflects all our class fields that we want to be settable via config file
	bool self_reflect(cgv::reflect::reflection_handler &srh)
	{
		return
			srh.reflect_member("mesh_filename", mesh_filename) &&
			srh.reflect_member("invent_missing_colors", invent_missing_colors) &&
			srh.reflect_member("show_surface", show_surface) &&
			srh.reflect_member("show_vertices", show_vertices) &&
			srh.reflect_member("show_wireframe", show_wireframe);
	}

	/// define all GUI elements for our mesh viewer
	void create_gui()
	{
		add_decorator("Mesh", "heading", "level=2");
		add_member_control(
			this, "invent missing per-vertex colors", invent_missing_colors, "check",
			"tooltip='After loading, invent per-vertex colors if the mesh did not include its own'"
		);
		add_gui(
			"mesh file ", mesh_filename, "file_name",
			"title='Load mesh from file';filter='mesh (obj):*.obj|all files:*.*';w=128"
		);
		connect_copy(
			add_button("Generate Conway polyhedron",
			           "tooltip='Replaces the current mesh with a procedural Conway polyhedron'")->click,
			cgv::signal::rebind(this, &mesh_viewer::create_conway_polyhedron)
		);

		add_decorator("", "separator");

		add_decorator("Display Settings", "heading", "level=2");
		bool show = begin_tree_node("vertices", show_vertices, false, "options='w=100';align=' '");
		add_member_control(this, "show", show_vertices, "toggle", "w=42;shortcut='w'", " ");
		add_member_control(this, "", sphere_style.surface_color, "", "w=42");
		if (show) {
			align("\a");
			add_gui("style", sphere_style);
			align("\b");
			end_tree_node(show_wireframe);
		}
		show = begin_tree_node("wireframe", show_wireframe, false, "options='w=100';align=' '");
		add_member_control(this, "show", show_wireframe, "toggle", "w=42;shortcut='w'", " ");
		add_member_control(this, "", cone_style.surface_color, "", "w=42");
		if (show) {
			align("\a");
			add_gui("style", cone_style);
			align("\b");
			end_tree_node(show_wireframe);
		}

		show = begin_tree_node("surface", show_surface, false, "options='w=100';align=' '");
		add_member_control(this, "show", show_surface, "toggle", "w=42;shortcut='s'", " ");
		add_member_control(this, "", surface_color, "", "w=42");
		if (show) {
			align("\a");
			add_member_control(this, "cull mode", cull_mode, "dropdown", "enums='none,back,front'");
			if (begin_tree_node("color_mapping", color_mapping)) {
				align("\a");
				add_gui("color mapping", color_mapping, "bit_field_control",
					"enums='COLOR_FRONT=1,COLOR_BACK=2,OPACITY_FRONT=4,OPACITY_BACK=8'");
				align("\b");
				end_tree_node(color_mapping);
			}
			add_member_control(this, "surface color", surface_color);
			add_member_control(this, "illumination", illumination_mode, "dropdown", "enums='none,one sided,two sided'");
			// this is how to add a ui for the materials read from an obj material file
			for (unsigned mi = 0; mi < mesh_info.ref_materials().size(); ++mi) {
				if (begin_tree_node(mesh_info.ref_materials()[mi]->get_name(), *mesh_info.ref_materials()[mi])) {
					align("\a");
					add_gui("mat", static_cast<cgv::media::illum::textured_surface_material&>(*mesh_info.ref_materials()[mi]));
					align("\b");
					end_tree_node(*mesh_info.ref_materials()[mi]);
				}
			}
			align("\b");
			end_tree_node(show_surface);
		}
		add_decorator("", "separator");

		add_decorator("Heightfield", "heading", "level=2");
		add_member_control(this, "visible view", visible_view, "value_slider", "min=0;max=44;ticks=true");
		add_member_control(this, "warping mode", warping_mode, "dropdown", "enums='simple interpolation, simple warping'");
	}
};

/// register the mesh_viewer drawable
#include <cgv/base/register.h>
cgv::base::object_registration<mesh_viewer> reg_mesh_viewer("");
cgv::base::registration_order_definition ro_def("holo_view_interactor;mesh_viewer");
