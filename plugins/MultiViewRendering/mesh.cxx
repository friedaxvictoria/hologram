#include <cgv/base/node.h>
#include <cgv/defines/quote.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/stereo_view.h>
#include <cgv/render/drawable.h>
#include <cgv/render/clipped_view.h>
#include <cgv_gl/gl/gl.h>
#include <cgv_gl/sphere_renderer.h>
#include <cgv_gl/cone_renderer.h>
#include <cgv_gl/gl/mesh_render_info.h>
#include <cgv_gl/box_wire_render_data.h>
#include <cgv_gl/gl/gl_context.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/application.h>
#include <cgv/media/color_scale.h>
#include <cgv/media/mesh/simple_mesh.h>
#include <cgv/gui/event_handler.h>
#include <cgv/gui/key_event.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/dialog.h>
#include <cgv/math/ftransform.h>
#include <cgv/render/render_buffer.h>
#include <cgv/render/managed_frame_buffer.h>
#include <cgv/render/clipped_view.h>
#include "tessellator.h"
#include "utilities.h"
#include "mesh.h"

using namespace cgv::base;
using namespace cgv::signal;
using namespace cgv::gui;
using namespace cgv::math;
using namespace cgv::render;
using namespace cgv::utils;
using namespace cgv::media::illum;

/*
This example illustrates how to use cgv::media::mesh::simple_mesh<T> and cgv::render::mesh_render_info
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

/// the constructor
mesh_viewer::mesh_viewer() : node("mesh_viewer")
{
	cull_mode = CM_BACKFACE;
	color_mapping = cgv::render::CM_COLOR;
	surface_color = rgb(0.75f, 0.25f, 1.0f);
	illumination_mode = IM_ONE_SIDED;

	// sphere_style.surface_color = rgb(0.8f, 0.3f, 0.3f);
	cone_style.surface_color = rgb(0.6f, 0.5f, 0.4f);

	eye_distance = 0.3f;

	nr_rendered_views_change();
}

/// reflect the name of our class
std::string mesh_viewer::get_type_name() const { return "mesh_viewer"; }

/// called when an instance of this class is registered with the Framework
void mesh_viewer::on_register()
{
	// --NOTE-- obviously, a window can only have one title, so it should be
	//  in the most central, "main"
	//          class that the project has (which this one, as it just provides some geometry to display, might
	//          not be anymore at some point)
	cgv::gui::application::get_window(0)->set("title", "Multi-View Rendering for Holographic Displays");
}

/// helper function that acts as a single point where processing of the mesh into a renderable form is happening
void mesh_viewer::process_mesh_for_rendering(context& ctx, bool update_view = true)
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
	M_bbox = M.compute_box();
	M_bbox_rd.clear();
	M_bbox_rd.add(M_bbox.get_center(), M_bbox.get_extent(), rgb(0.75f));

	// [re-]compute mesh render info
	mesh_info.destruct(ctx);
	mesh_info.construct(ctx, M);
	// bind mesh attributes to standard surface shader program
	mesh_info.bind(ctx, ctx.ref_surface_shader_program(true), true);
	mesh_info.bind_wireframe(ctx, ref_cone_renderer(ctx).ref_prog(), true);

	// [re-]compute mesh render info for the geometry shader acceleration
	mesh_for_geo_info.destruct(ctx);
	mesh_for_geo_info.construct(ctx, M);
	// bind mesh attributes to standard surface shader program
	mesh_for_geo_info.bind(ctx, test.geometry_shader, true);
	mesh_for_geo_info.bind_wireframe(ctx, ref_cone_renderer(ctx).ref_prog(), true);

	// [re-]compute mesh render info for visualising the holes
	mesh_for_holes_info.destruct(ctx);
	mesh_for_holes_info.construct(ctx, M);
	mesh_for_holes_info.bind(ctx, test.holes_shader, true);

	// update sphere attribute array manager
	sphere_renderer& sr = ref_sphere_renderer(ctx);
	sr.enable_attribute_array_manager(ctx, sphere_aam);
	sr.set_position_array(ctx, M.get_positions());
	if (M.has_colors())
		sr.set_color_array(ctx, *reinterpret_cast<const std::vector<rgb>*>(M.get_color_data_vector_ptr()));

	// adjust camera parameters when requested
	update_view_after_mesh_processed |= update_view;

	// update view
	// make sure we have the view available
	if (!view) {
		view = dynamic_cast<stereo_view*>(find_view_as_node());
		// dynamic_cast<node*>(view)->set(
		//"clip_relative_to_extent",
		//   true); // comment this line to use default behaviour of adaptating znear/zfar to scene
	}
	// focus view on new mesh
	view->set_scene_extent(M_bbox);
	view->set_focus(M_bbox.get_center());
	view->set_y_extent_at_focus(M_bbox.get_extent().length());

	// ensure that materials are presented in gui
	post_recreate_gui();
}

/// helper function that will make sure we have a per-vertex color attribute
void mesh_viewer::ensure_mesh_colors()
{
	if (M.has_colors())
		return;
	const int nr_positions = M.get_nr_positions();
	M.ensure_colors(cgv::media::CT_RGB, nr_positions);
	double dummy;
#pragma omp for
	for (int i = 0; i < nr_positions; i++) {
		M.set_color(i, cgv::media::color_scale(modf(double(20 * i) / double(nr_positions - 1), &dummy)));
	}
}

/// react to mouse & keyboard events
bool mesh_viewer::handle(cgv::gui::event& e)
{
	if (e.get_kind() == cgv::gui::EID_KEY) {
		const auto& ke = static_cast<cgv::gui::key_event&>(e);
		const auto ka = ke.get_action();

		switch (ke.get_key()) {
		case Keys::KEY_Enter:
			if (ka == KeyAction::KA_PRESS && !test.shoot_heightmap) {
				test.shoot_heightmap = true;
				post_redraw();
				return true;
			}
		default:
			  /* DoNothing() */;
		}
	}

	return false;
}

/// output help text for keyboard shortcuts
void mesh_viewer::stream_help(std::ostream& os)
{
	os << "mesh_viewer:" << std::endl
	   << "\ttoggle mesh surface[s], toggle mesh vertices[v], toggle mesh wireframe[w]," << std::endl
	   << "\ttoggle mesh bounding box[b], toggle heightmap[h], capture cur. heightmap[ENTER]" << std::endl
	   << "\ttoggle heightmap pruning[p]" << std::endl;
}

/// react to our class fields being set via the GUI or via reflection (e.g. from a config file)
void mesh_viewer::on_set(void* member_ptr)
{
	if (member_ptr == &mesh_filename) {
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

	////
	// SECTION: 3D image warping baseline test

	// react to changed heightmap oversampling factor
	if (member_ptr == &test.heightmap_oversampling)
		test.render_fbo->set_size({-1, -1}); // causes FBO and heightmap regeneration during ::init_frame()

	// END: 3D image warping baseline test
	////

	update_member(member_ptr);
	post_redraw();
}

/// clears whatever mesh is currently loaded an creates a Conway polyhedron instead
void mesh_viewer::create_conway_polyhedron()
{
	M.clear();
	M.construct_conway_polyhedron("dtI");
	mesh_filename = "<generated>";
	process_mesh_for_rendering(*get_context());
}

/// perform initialization that can not be done in the constructor as it requires a fully functional
/// graphics context
bool mesh_viewer::init(context& ctx)
{
	// init render components
	ref_sphere_renderer(ctx, 1);
	ref_cone_renderer(ctx, 1);
	bool success = sphere_aam.init(ctx);
	success &= cone_aam.init(ctx);
	success &= M_bbox_rd.init(ctx);

	// in the blank state (without anything loaded), we just display a simple Conway polyhedron
	create_conway_polyhedron();

	////
	// SECTION: 3D image warping baseline test

	// create our offscreem framebuffers used to test image warping
	for (unsigned i = 0; i < 3; i++) {
		test.render_fbo[i].add_attachment("depth", "[D]");
		test.render_fbo[i].add_attachment("color", "uint8[R,G,B,A]");
	}

	// init our baseline method test shaders
	success &= test.baseline_shader.build_program(ctx, "baseline.glpr", true);
	success &= test.warp_shader.build_program(ctx, "warp.glpr", true);
	success &= test.holes_shader.build_program(ctx, "holes.glpr", true);
	success &= test.geometry_shader.build_program(ctx, "geometry.glpr", true);

	test.geometry_shader.specify_standard_uniforms(true, true, true, true);
	test.geometry_shader.specify_standard_vertex_attribute_names(ctx, true, true, true);
	test.geometry_shader.allow_context_to_set_color(true);

	// END: 3D image warping baseline test
	////

	// report success (or lack thereof)
	return success;
}

/// unload any currently loaded mesh data
void mesh_viewer::clear(context& ctx)
{
	M_bbox_rd.destruct(ctx);
	ref_cone_renderer(ctx, -1);
	ref_sphere_renderer(ctx, -1);
	sphere_aam.destruct(ctx);
	cone_aam.destruct(ctx);
}

/// perform any kind of operation that should take place before all drawables start executing their ::draw()
/// methods
void mesh_viewer::init_frame(context& ctx)
{
	////
	// SECTION: 3D image warping baseline test

	// check if render fbo needs re-initialization (we always change all 3 at the same time, so we just check the
	// center one)
	auto& fb = test.render_fbo[1].ref_frame_buffer();
	if (fb.get_width() != (int)float(ctx.get_width() * test.heightmap_oversampling) ||
		fb.get_height() != (int)float(ctx.get_height() * test.heightmap_oversampling))
	{
		// set resolution according to oversampling factor
		test.render_fbo[0].set_size({(int)float(ctx.get_width() * test.heightmap_oversampling),
									 (int)float(ctx.get_height() * test.heightmap_oversampling)});
		test.render_fbo[0].ensure(ctx);
		test.render_fbo[1].set_size(test.render_fbo[0].get_size());
		test.render_fbo[1].ensure(ctx);
		test.render_fbo[2].set_size(test.render_fbo[1].get_size());
		test.render_fbo[2].ensure(ctx);

		// (re-)tessellate our heightmap
		// - get view frustum information
		const auto res = test.render_fbo[1].get_size();
		const float half_aspect = (float)res.x() / (2 * res.y());
		// - tesselate
		test.heightmap_baseline = tessellator::quad(ctx, test.baseline_shader, {-half_aspect, -.5f, .0f},
													{half_aspect, .5f, .0f}, res.x(), res.y());
		test.heightmap_warp = tessellator::quad(ctx, test.warp_shader, {-half_aspect, -.5f, .0f},
												{half_aspect, .5f, .0f}, res.x(), res.y());

		// make sure we re-shoot the heightmap
		test.shoot_heightmap = true;
	}

	// END: 3D image warping baseline test
	////
}

void mesh_viewer::draw_holes(context& ctx)
{
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);

	mesh_for_holes_info.draw_all(ctx);

	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
}

void mesh_viewer::set_params_for_gemoetry(float zero_parallax, float eye_distance, float eye, float num_holo_views)
{
	test.zero_parallax = zero_parallax;
	test.eye_distance = eye_distance;
	test.eye = eye;
	test.num_holo_views = num_holo_views;
	test.with_geometry = true;
}

void mesh_viewer::draw_geometry_shader(context& ctx)
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

	test.geometry_shader.set_uniform(ctx, "eye_sep", test.eye_distance);
	test.geometry_shader.set_uniform(ctx, "eye", test.eye);
	test.geometry_shader.set_uniform(ctx, "num_holo_views", test.num_holo_views);
	test.geometry_shader.set_uniform(ctx, "zero_parallax", test.zero_parallax);
	test.geometry_shader.set_uniform(ctx, "culling_mode", (int)cull_mode);
	test.geometry_shader.set_uniform(ctx, "map_color_to_material", (int)color_mapping);
	test.geometry_shader.set_uniform(ctx, "illumination_mode", (int)illumination_mode);

	mesh_for_geo_info.draw_all(ctx);

	// recover opengl culling mode
	if (is_culling)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);
	glCullFace(cull_face);

	test.with_geometry = false;
}

/// draw the mesh surface
void mesh_viewer::draw_surface(context& ctx, bool opaque_part)
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
	shader_program& prog = ctx.ref_surface_shader_program(true);
	prog.set_uniform(ctx, "culling_mode", (int)cull_mode);
	prog.set_uniform(ctx, "map_color_to_material", (int)color_mapping);
	prog.set_uniform(ctx, "illumination_mode", (int)illumination_mode);
	// set default surface color for color mapping which only affects
	// rendering if mesh does not have per vertex colors and color_mapping is on
	// prog.set_attribute(ctx, prog.get_color_index(), surface_color);
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
void mesh_viewer::draw(context& ctx)
{
	////
	// SECTION: 3D image warping baseline test

	// --NOTE-- bind appropriate render_fbo here (e.g. just render_fbo[1] for the
	//          initial mono depth map test case)
	/* test.render_heightmap = true;

	if (test.shoot_heightmap) {
		static const vec3 right_local(1, 0, 0);
		const mat4 invMV = inv(ctx.get_modelview_matrix());
		cam_dir = view->get_view_dir();
		const vec3 cam_pos = view->get_eye(), cam_right = normalize(w_clip(invMV.mul_pos(right_local)) - cam_pos),
	cam_up = view->get_view_up_dir(); // --CAREFUL-- get_view_up_dir is not 100% reliable,
													 // consider inferring from modelview matrix!
		// - find out where the scene ends along world-space camera viewing direction
		auto [bmin, bmax] = project_box_onto_dir(M_bbox, cam_dir);
		// - find out how large the quad needs to be to fill the whole frustum at that point
		const float cam_depth_world = dot(cam_pos, cam_dir), heightmap_depth_eye = bmax - cam_depth_world,
					aspect = (float)ctx.get_width() / ctx.get_height();
		test.y_ext = (float)view->get_y_extent_at_depth(heightmap_depth_eye, false);
		test.x_ext = test.y_ext * aspect;
		// - position the heightmap just behind the scene, facing the camera at the moment of capture:
		//   M = Transl(cam_pos + cam_dir*heightmap_depth_eye) * Rot(cam_right, cam_up, cam_dir) * Scale(y_ext)
		test.heightmap_trans = mat4({vec4(test.y_ext * cam_right, 0), vec4(test.y_ext * cam_up, 0), vec4(test.y_ext *
	cam_dir, 0), vec4(cam_pos + cam_dir * heightmap_depth_eye, 1)}); test.modelview_source = ctx.get_modelview_matrix()
	* test.heightmap_trans;

		const float x_ext_half = 0.5f * test.x_ext, y_ext_half = 0.5f * test.y_ext;
		test.znear = bmin - cam_depth_world;

		// change projection matrix to orthogonal according to the current (view-dependent) extends of our
		// heightmap or use perspective projection matrix
		if (!test.ortho)
			test.proj_for_render = ctx.get_projection_matrix();
		else
			test.proj_for_render = ortho4(-x_ext_half, x_ext_half, -y_ext_half, y_ext_half, test.znear,
	heightmap_depth_eye);
	}

	for (int i = 0; i < test.nr_renders; i++) {
		if (test.shoot_heightmap || (test.nr_rendered_views_changed && test.render_heightmap)) {
			test.render_fbo[i].enable(ctx);
			// make the heightmap quad very slightly visible even where it doesn't contain scene geometry
			// glClearColor(.03125f, .03125f, .03125f, 1.f);

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			// move the heightmap to the correct place in world-space
			// - determine the camera orientation

			float x_of_cam = 0.5f * test.render_offset[i] * view->get_eye_distance() * test.x_ext;
			test.eye_source[i] = vec3(x_of_cam, 0, 0);

			mat4 shear, translate;
			shear.identity();
			shear(0, 2) = x_of_cam / view->get_depth_of_focus();

			translate.identity();
			translate(0, 3) = x_of_cam;

			test.inv_mat_proj_render[i] = inv(test.proj_for_render * shear * translate);
			std::cout << "proj in mesh: " << i << ", "<< test.inv_mat_proj_render[i] << std::endl;

			// test.p_1[i] = compute_frustum_model_image_warp(-test.x_ext/2, test.x_ext/2, -test.y_ext/2,
			// test.y_ext, ctx.get_width(), ctx.get_height(), test.znear);

			ctx.push_projection_matrix();
			ctx.set_projection_matrix(inv(test.inv_mat_proj_render[i]));
			ctx.push_modelview_matrix();
			ctx.set_modelview_matrix(test.modelview_source * inv(test.heightmap_trans));
		}



		// END: 3D image warping baseline test
		////

		 if (show_vertices && (test.shoot_heightmap || !test.render_heightmap))
		{
			sphere_renderer& sr = ref_sphere_renderer(ctx);
			sr.set_render_style(sphere_style);
			sr.enable_attribute_array_manager(ctx, sphere_aam);
			sr.render(ctx, 0, M.get_nr_positions());
			sr.disable_attribute_array_manager(ctx, sphere_aam);
		}
		if (show_wireframe && (test.shoot_heightmap || !test.render_heightmap)) {
			cone_renderer& cr = ref_cone_renderer(ctx);
			cr.set_render_style(cone_style);
			if (cr.enable(ctx)) {
				mesh_info.draw_wireframe(ctx);
				cr.disable(ctx);
			}
		}
		if (show_surface && (test.shoot_heightmap || !test.render_heightmap || test.nr_rendered_views_changed))
			draw_surface(ctx, false); // --NOTE-- set parameter to true once done testing the image warp to
									  // re-enable correct handling of transparent mesh parts

		// draw the mesh bounding box if we're not currently capturing the heightmap
		if (show_bbox && !test.shoot_heightmap)
			M_bbox_rd.render(ctx, ref_box_wire_renderer(ctx), box_wire_render_style());

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

		// we only position the heightmap in space when it's being "shot", it will remain there afterwards
		// so it can be inspected from all sides and compared to the "real" scene
		if (test.shoot_heightmap || (test.nr_rendered_views_changed && test.render_heightmap)) {
			// disable the offscreen framebuffer and reset projection matrix
			test.render_fbo[i].disable(ctx);
			ctx.pop_projection_matrix();
			ctx.pop_modelview_matrix();

			// shoot complete
			if (i == test.nr_renders-1) {
				test.shoot_heightmap = false;
				test.nr_rendered_views_changed = false;
			}

			// post a redraw so the scene is rendered once more into the main framebuffer (at this point it will
			// have only been drawn into the offscreen framebuffer)
			post_redraw();
		}

	}

	if (test.render_heightmap) {

		test.eye_target =
			  vec3(0.5f * (2.0f*test.visible_view / 44.0 - 1.0f) * view->get_eye_distance() * test.x_ext, 0, 0);

		mat3 p_2 = compute_frustum_model_image_warp(-test.x_ext / 2, test.x_ext / 2, -test.y_ext / 2, test.y_ext,
													ctx.get_width(), ctx.get_height(), test.znear);

		if (test.show_holes) {
			// render pass for displaying holes when warping in green
			glDisable(GL_CULL_FACE);
			glDisable(GL_DEPTH_TEST);

			mesh_for_holes_info.draw_all(ctx);

			glEnable(GL_CULL_FACE);
			glEnable(GL_DEPTH_TEST);
		}

		if (test.warp_mode == test.BASELINE) {
			for (int i = 0; i < test.nr_renders; i++) {

				texture &color_tex = *test.render_fbo[i].attachment_texture_ptr("color"),
						&depth_tex = *test.render_fbo[i].attachment_texture_ptr("depth");
				// render pass for baseline approach
				color_tex.enable(ctx, 0);
				test.baseline_shader.set_uniform(ctx, "color", 0);
				depth_tex.enable(ctx, 1);
				test.baseline_shader.set_uniform(ctx, "depth", 1);

				test.baseline_shader.set_uniform(ctx, "inv_proj_source", test.inv_mat_proj_render[i]);
				test.baseline_shader.set_uniform(ctx, "modelview_source", test.modelview_source);
				test.baseline_shader.set_uniform(ctx, "prune_empty", test.prune_heightmap);
				test.baseline_shader.set_uniform(ctx, "with_interpolated_holes",
												 test.with_interpolated_holes && test.render_offset[i] == 0.0f);
				test.baseline_shader.set_uniform(ctx, "epsilon", test.epsilon);

				test.baseline_shader.enable(ctx);
				glDisable(GL_CULL_FACE);
				ctx.push_modelview_matrix();
				ctx.mul_modelview_matrix(test.heightmap_trans);
				test.heightmap_baseline.draw(ctx);
				ctx.pop_modelview_matrix();
				glEnable(GL_CULL_FACE);
				test.baseline_shader.disable(ctx);
				depth_tex.disable(ctx);
				color_tex.disable(ctx);
			}
		}
		else {
			for (int i = 0; i < test.nr_renders; i++) {

				texture &color_tex = *test.render_fbo[i].attachment_texture_ptr("color"),
						&depth_tex = *test.render_fbo[i].attachment_texture_ptr("depth");

				color_tex.enable(ctx, 0);
				test.warp_shader.set_uniform(ctx, "color", 0);
				depth_tex.enable(ctx, 1);
				test.warp_shader.set_uniform(ctx, "depth", 1);

				test.warp_shader.set_uniform(ctx, "p_1", test.p_1[i]);
				test.warp_shader.set_uniform(ctx, "p_2", p_2);
				test.warp_shader.set_uniform(ctx, "inv_mvp_source", inv(test.modelview_source) *
	test.inv_mat_proj_render[i]); test.warp_shader.set_uniform(ctx, "eye_source", test.eye_source[i]);
				test.warp_shader.set_uniform(ctx, "eye_target", test.eye_target);
				test.warp_shader.set_uniform(ctx, "width", (float)(ctx.get_width() * test.heightmap_oversampling));
				test.warp_shader.set_uniform(ctx, "height", (float)(ctx.get_height() * test.heightmap_oversampling));
				test.warp_shader.set_uniform(ctx, "fac",
											 fac);

				test.warp_shader.enable(ctx);
				glDisable(GL_CULL_FACE);
				ctx.push_modelview_matrix();
				ctx.mul_modelview_matrix(test.heightmap_trans);
				test.heightmap_warp.draw(ctx);
				ctx.pop_modelview_matrix();
				glEnable(GL_CULL_FACE);
				test.warp_shader.disable(ctx);
				depth_tex.disable(ctx);
				color_tex.disable(ctx);
			}
		}
	} */

	if (show_vertices) {
		sphere_renderer& sr = ref_sphere_renderer(ctx);
		sr.set_render_style(sphere_style);
		sr.enable_attribute_array_manager(ctx, sphere_aam);
		sr.render(ctx, 0, M.get_nr_positions());
		sr.disable_attribute_array_manager(ctx, sphere_aam);
	}
	if (show_wireframe) {
		cone_renderer& cr = ref_cone_renderer(ctx);
		cr.set_render_style(cone_style);
		if (cr.enable(ctx)) {
			mesh_info.draw_wireframe(ctx);
			cr.disable(ctx);
		}
	}
	if (show_surface)
		if (test.with_geometry)
			draw_geometry_shader(ctx);
		else
			draw_surface(ctx, true); // --NOTE-- set parameter to true once done testing the image warp to
								 // re-enable correct handling of transparent mesh parts

	// draw the mesh bounding box if we're not currently capturing the heightmap
	if (show_bbox)
		M_bbox_rd.render(ctx, ref_box_wire_renderer(ctx), box_wire_render_style());

	// END: 3D image warping baseline test
	////
}

/// perform any kind of operation that should take place after all drawables have executed their ::draw()
/// methods
void mesh_viewer::finish_frame(context& ctx)
{
	if (show_surface)
		draw_surface(ctx, false); // --NOTE-- uncomment once done testing the image warp to re-enable correct handling
								  // of transparent mesh parts
}

/// reflects all our class fields that we want to be settable via config file
bool mesh_viewer::self_reflect(cgv::reflect::reflection_handler& srh)
{
	return srh.reflect_member("mesh_filename", mesh_filename) &&
		   srh.reflect_member("invent_missing_colors", invent_missing_colors) &&
		   srh.reflect_member("show_bbox", show_bbox) && srh.reflect_member("show_surface", show_surface) &&
		   srh.reflect_member("show_vertices", show_vertices) && srh.reflect_member("show_wireframe", show_wireframe) &&
		   srh.reflect_member("test__render_heightmap", test.render_heightmap) &&
		   srh.reflect_member("test__warp_mode", test.warp_mode) &&
		   srh.reflect_member("test__heightmap_oversampling", test.heightmap_oversampling) &&
		   srh.reflect_member("test__epsilon", test.epsilon) &&
		   srh.reflect_member("test__nr_renders", test.nr_renders) &&
		   srh.reflect_member("test__prune_heightmap", test.prune_heightmap) &&
		   srh.reflect_member("test__with_interpolated_holes", test.with_interpolated_holes) &&
		   srh.reflect_member("test__show_holes", test.show_holes) && srh.reflect_member("test__ortho", test.ortho);
}

/// define all GUI elements for our mesh viewer
// - helper one-shot methods needed as button callbacks (until the framework properly supports C++ lambdas)
void mesh_viewer::on_shoot(void)
{
	test.shoot_heightmap = true;
	post_redraw();
}

void mesh_viewer::nr_rendered_views_change(void)
{
	test.nr_rendered_views_changed = true;
	test.render_offset[1] = -1.0;
	test.render_offset[2] = 1.0;
	if (test.nr_renders == 1 || test.nr_renders == 3) {
		test.render_offset[0] = 0.0;
	}
	else {
		test.render_offset[0] = 1.0;
	}
	post_redraw();
}

// - the actual method
void mesh_viewer::create_gui()
{
	add_decorator("Mesh", "heading", "level=2");
	add_member_control(this, "invent missing per-vertex colors", invent_missing_colors, "check",
					   "tooltip='After loading, invent per-vertex colors if the mesh did not include its own'");
	add_gui("mesh file ", mesh_filename, "file_name",
			"title='Load mesh from file';filter='mesh (obj):*.obj|all files:*.*';w=128");
	connect_copy(add_button("Generate Conway polyhedron",
							"tooltip='Replaces the current mesh with a procedural Conway polyhedron'")
					   ->click,
				 cgv::signal::rebind(this, &mesh_viewer::create_conway_polyhedron));

	add_decorator("", "separator");

	add_decorator("Display Settings", "heading", "level=2");
	bool show = begin_tree_node("vertices", show_vertices, false, "options='w=100';align=' '");
	add_member_control(this, "show", show_vertices, "toggle", "w=42;shortcut='v'", " ");
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
				add_gui("mat",
						static_cast<cgv::media::illum::textured_surface_material&>(*mesh_info.ref_materials()[mi]));
				align("\b");
				end_tree_node(*mesh_info.ref_materials()[mi]);
			}
		}
		align("\b");
		end_tree_node(show_surface);
	}
	add_member_control(this, "show bounding box", show_bbox, "check", "shortcut='b'");

	////
	// SECTION: 3D image warping baseline test

	add_decorator("", "separator");
	add_decorator("Stereo Viewing", "heading", "level=2");

	add_decorator("", "separator");
	add_decorator("Baseline 3D Image Warp", "heading", "level=2");
	add_member_control(this, "render as heightmap", test.render_heightmap, "toggle", "shortcut='h'");
	add_member_control(this, "Warp Mode", test.warp_mode, "dropdown", "enums='baseline, image_warp'");
	add_member_control(this, "oversampling factor", test.heightmap_oversampling, "value_slider",
					   "min=0.5;max=4;step=0.5;tooltip='multiply viewport resolution by this factor when "
					   "determining heightmap resolution'");
	add_member_control(this, "fac", fac, "value_slider",
					   "min=0;max=0.5;step=0.00001;tooltip='how strong should the discarding of cliff artefacts be'");
	add_member_control(this, "epsilon", test.epsilon, "value_slider",
					   "min=0;max=0.1;step=0.00001;tooltip='how strong should the discarding of cliff artefacts be'");
	add_member_control(this, "visible view", test.visible_view, "value_slider",
					   "min=0;max=44;step=1;tooltip='change how many views are rendered properly'");
	connect_copy(add_member_control(this, "nr rendered views", test.nr_renders, "value_slider",
									"min=1;max=3;step=1;tooltip='change how many views are rendered properly'")
					   ->value_change,
				 cgv::signal::rebind(this, &mesh_viewer::nr_rendered_views_change));
	add_member_control(this, "prune empty areas", test.prune_heightmap, "check",
					   "tooltip='discard heightmap fragments which dont represent valid geometry';shortcut='p'");
	add_member_control(this, "show holes", test.show_holes, "check");
	add_member_control(this, "smooth holes into background", test.with_interpolated_holes, "check");
	connect_copy(add_member_control(this, "with orthographic matix", test.ortho, "check")->value_change,
				 cgv::signal::rebind(this, &mesh_viewer::on_shoot));
	connect_copy(add_button("Shoot!", "tooltip='Updates the heightmap from the current view'")->click,
				 cgv::signal::rebind(this, &mesh_viewer::on_shoot));

	// END: 3D image warping baseline test
	////
}

/// register the mesh_viewer drawable
#include <cgv/base/register.h>
cgv::base::object_registration<mesh_viewer> reg_mesh_viewer("");
cgv::base::registration_order_definition ro_def("holo_view_interactor;mesh_viewer");