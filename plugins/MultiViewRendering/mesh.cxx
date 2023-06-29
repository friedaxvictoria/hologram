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

	sphere_style.surface_color = rgb(0.8f, 0.3f, 0.3f);
	cone_style.surface_color = rgb(0.6f, 0.5f, 0.4f);

	eye_distance = 0.3f;
}

/// reflect the name of our class
std::string mesh_viewer::get_type_name() const { return "mesh_viewer"; }

/// called when an instance of this class is registered with the Framework
void mesh_viewer::on_register()
{
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
	mesh_for_geo_info.bind(ctx, geometry_shader, true);
	mesh_for_geo_info.bind_wireframe(ctx, ref_cone_renderer(ctx).ref_prog(), true);

	// [re-]compute mesh render info for visualising the holes
	mesh_for_holes_info.destruct(ctx);
	mesh_for_holes_info.construct(ctx, M);
	mesh_for_holes_info.bind(ctx, holes_shader, true);

	// update sphere attribute array manager
	sphere_renderer& sr = ref_sphere_renderer(ctx);
	sr.enable_attribute_array_manager(ctx, sphere_aam);
	sr.set_position_array(ctx, M.get_positions());
	if (M.has_colors())
		sr.set_color_array(ctx, *reinterpret_cast<const std::vector<rgb>*>(M.get_color_data_vector_ptr()));
	else {
		std::vector<rgb> colour;
		colour.resize(M.get_positions().size());
		for (size_t i = 0; i < M.get_positions().size(); ++i) {
			colour[i] = rgb(sphere_style.surface_color);
		}
		sr.set_color_array(ctx, colour);
	}

	// adjust camera parameters when requested
	update_view_after_mesh_processed |= update_view;

	// update view
	// make sure we have the view available
	if (!view)
		view = dynamic_cast<stereo_view*>(find_view_as_node());
	focus_mesh();

	// ensure that materials are presented in gui
	post_recreate_gui();
}

void mesh_viewer::focus_mesh() {
	// focus view on new mesh
	view->set_scene_extent(M_bbox);
	view->set_focus(M_bbox.get_center());
	view->set_y_extent_at_focus(M_bbox.get_extent().length());
}

int mesh_viewer::get_number_positions() { 
	return M.get_nr_positions(); 
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
		double v = modf(double(20 * i) / double(nr_positions - 1), &dummy);
		// interpolate between ##0072BD, #7E2F8E, and #4DBEEE
		if (v<0.5) M.set_color(
			  i, cgv::media::color<float, cgv::media::RGB>(0.4940 * 2 * v, 0.4470 * (1 - 2 * v) + 0.1840 * 2 * v,
																  0.7410 * (1 - 2 * v) + 0.5560 * 2 * v));
		else M.set_color(i, cgv::media::color<float, cgv::media::RGB>(0.4940 * (2 - 2 * v) + 0.3010 * (2 * v - 1),
																	   0.1840 * (2 - 2 * v) + 0.7450 * (2 * v - 1),
																	 0.5560 * (2 - 2 * v) + 0.9330 * (2 * v - 1)));
	}
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
	if (member_ptr == &sphere_style.surface_color && !M.has_colors())
		process_mesh_for_rendering(*get_context());

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

	success &= holes_shader.build_program(ctx, "holes.glpr", true);
	success &= geometry_shader.build_program(ctx, "geometry.glpr", true);

	geometry_shader.specify_standard_uniforms(true, true, true, true);
	geometry_shader.specify_standard_vertex_attribute_names(ctx, true, true, true);
	geometry_shader.allow_context_to_set_color(true);

	// report success (or lack thereof)
	return success;
}

bool mesh_viewer::handle(cgv::gui::event& e) { return true; };

/// unload any currently loaded mesh data
void mesh_viewer::clear(context& ctx)
{
	M_bbox_rd.destruct(ctx);
	ref_cone_renderer(ctx, -1);
	ref_sphere_renderer(ctx, -1);
	sphere_aam.destruct(ctx);
	cone_aam.destruct(ctx);
}

void mesh_viewer::init_frame(context& ctx){}

// draw holes in yellow behind the rendered object
void mesh_viewer::draw_holes(context& ctx)
{
	glDisable(GL_CULL_FACE);
	//disable depth test so that hole visualisation is always in the back
	glDisable(GL_DEPTH_TEST);

	mesh_for_holes_info.draw_all(ctx);

	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
}

void mesh_viewer::set_params_for_gemoetry(float zero_parallax, float eye_distance, float eye, float num_holo_views)
{
	mesh_viewer::zero_parallax = zero_parallax;
	mesh_viewer::eye_distance = eye_distance;
	mesh_viewer::eye = eye;
	mesh_viewer::num_holo_views = num_holo_views;
	with_geometry = true;
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

	geometry_shader.set_uniform(ctx, "eye_sep", eye_distance);
	geometry_shader.set_uniform(ctx, "eye", eye);
	geometry_shader.set_uniform(ctx, "view_offset", 2/float(num_holo_views));
	geometry_shader.set_uniform(ctx, "zero_parallax", zero_parallax);
	geometry_shader.set_uniform(ctx, "culling_mode", (int)cull_mode);
	geometry_shader.set_uniform(ctx, "map_color_to_material", (int)color_mapping);
	geometry_shader.set_uniform(ctx, "illumination_mode", (int)illumination_mode);

	ctx.set_color(surface_color);

	mesh_for_geo_info.draw_all(ctx);

	// recover opengl culling mode
	if (is_culling)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);
	glCullFace(cull_face);
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

	// render with standard shader program or geometry shader program depending on what was selected in the gui
	if (show_surface) {
		if (with_geometry)
			draw_geometry_shader(ctx);
		else
			draw_surface(ctx, true);
	}

	if (show_bbox)
		M_bbox_rd.render(ctx, ref_box_wire_renderer(ctx), box_wire_render_style());
}

/// perform any kind of operation that should take place after all drawables have executed their ::draw()
/// methods
void mesh_viewer::finish_frame(context& ctx)
{
	if (show_surface)
		// render with standard shader program or geometry shader program depending on what was selected in the gui
		if (with_geometry) {
			draw_geometry_shader(ctx);
			with_geometry = false;

		}
		else
			draw_surface(ctx, false); 
}

/// reflects all our class fields that we want to be settable via config file
bool mesh_viewer::self_reflect(cgv::reflect::reflection_handler& srh)
{
	return srh.reflect_member("mesh_filename", mesh_filename) &&
		   srh.reflect_member("invent_missing_colors", invent_missing_colors) &&
		   srh.reflect_member("show_bbox", show_bbox) && srh.reflect_member("show_surface", show_surface) &&
		   srh.reflect_member("show_vertices", show_vertices) && srh.reflect_member("show_wireframe", show_wireframe);
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
		end_tree_node(show_vertices);
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
}

/// register the mesh_viewer drawable
#include <cgv/base/register.h>
cgv::base::object_registration<mesh_viewer> reg_mesh_viewer("");
cgv::base::registration_order_definition ro_def("holo_view_interactor;mesh_viewer");