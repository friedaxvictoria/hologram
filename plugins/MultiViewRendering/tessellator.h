#pragma once

#include <cgv/math/fvec.h>
#include <cgv/render/context.h>
#include <cgv/render/render_types.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/vertex_buffer.h>
#include <cgv/render/attribute_array_binding.h>
#include "lib_begin.h"


////
// utilities for tessellating a quad as a triangle mesh

/// struct encapsulating ready-to-render geometry stored on the GPU
struct CGV_API GPUgeometry
{
	cgv::render::vertex_buffer vbo;
	cgv::render::vertex_buffer indices;
	cgv::render::attribute_array_binding vao;
	unsigned restart_index; // non-zero if the geometry is a triangle strip, zero otherwise

	GPUgeometry() : indices(cgv::render::VBT_INDICES) {}
	GPUgeometry(const GPUgeometry &) = delete;
	GPUgeometry(GPUgeometry &&other)
		: vbo(std::move(other.vbo)), indices(std::move(other.indices)), vao(std::move(other.vao)),
		  restart_index(other.restart_index)
	{
		// prevent the moved-from render components from destroying the now-moved resources
		other.vbo.ctx_ptr = nullptr;
		other.indices.ctx_ptr = nullptr;
		other.vao.ctx_ptr = nullptr;
	}
	GPUgeometry& operator= (GPUgeometry&& other)
	{
		vbo = std::move(other.vbo);         other.vbo.ctx_ptr = nullptr;
		indices = std::move(other.indices); other.indices.ctx_ptr = nullptr;
		vao = std::move(other.vao);         other.vao.ctx_ptr = nullptr;
		restart_index = other.restart_index;
		return *this;
	}

	void draw (cgv::render::context &ctx);
};


/// static class providing methods to create ready-to-render geometry for simple shapes. For now, all
/// geometry defines positions and a 2D texture coordinates for the vertices.
class CGV_API tessellator
{
public:
	typedef cgv::render::render_types::vec2 vec2;
	typedef cgv::render::render_types::vec3 vec3;
	typedef cgv::render::render_types::vec4 vec4;

	/// tessellates a quad as a grid of num_samples0 � num_samples1 vertices connected in a triangle
	/// strip, and uploads the resulting geometry to the GPU associated with the given context, binding
	/// the vertex attributes to the "position" and "texcoord" inputs of the given shader program.
	/// dimension '0' refers to the larger extend of the quad, while '1' refers to the smaller one.
	static GPUgeometry quad (
		cgv::render::context &ctx, cgv::render::shader_program &shader, const vec3 &min_corner,
		const vec3 &max_corner, unsigned num_samples0, unsigned num_samples1
	);
};
