
#include <libs/cgv_gl/gl/gl.h>
#include <cgv/math/functions.h>
#include "tessellator.h"


using namespace cgv::math;
using namespace cgv::render;


// anonymous namespace for helper functions
namespace {
	// given a "mask" number m∈{0,1,2}, returns the other two (i.e. unmasked) numbers
	inline const std::pair<const unsigned, const unsigned> select_unmasked (const unsigned mask) {
		return {mask != 0 ? 0 : (mask != 1 ? 1 : 2),
		        mask != 2 ? 2 : 1};
	}

	// given a number i∈{0,1,2} "masking" a dimension from the given 3-vector, return the indices
	// of the other dimensions ordered according to the value of their corresponding components
	// (index of larger component first, index of smaller component second)
	template <class T>
	inline const std::pair<const unsigned, const unsigned> ordered_unmasked_dims (
		const cgv::math::fvec<T, 3> &vec, const unsigned mask
	) {
		auto [dim0, dim1] = select_unmasked(mask);
		if (vec(dim1) > vec(dim0))
			return {dim1, dim0};
		else
			return {dim0, dim1};
	}

	// the format of our vertices
	struct vertex {
		typedef typename tessellator::vec2 vec2;
		typedef typename tessellator::vec4 vec4;

		static const cgv::render::type_descriptor vec2desc;
		static const cgv::render::type_descriptor vec4desc;

		vec4 position;
		vec2 texcoord;

		static constexpr size_t position_offset = 0,
		                        texcoord_offset = sizeof(vec4);
	};
	const cgv::render::type_descriptor vertex::vec2desc = element_descriptor_traits<vertex::vec2>::get_type_descriptor(vertex::vec2());
	const cgv::render::type_descriptor vertex::vec4desc = element_descriptor_traits<vertex::vec4>::get_type_descriptor(vertex::vec4());
};


void GPUgeometry::draw(context& ctx)
{
	vao.enable(ctx);
	const GLsizei num_indices = (GLsizei)indices.get_size_in_bytes() / sizeof(unsigned);
	if (restart_index) {
		glEnable(GL_PRIMITIVE_RESTART);
		glPrimitiveRestartIndex(restart_index);
		glDrawElements(GL_TRIANGLE_STRIP, num_indices, GL_UNSIGNED_INT, nullptr);
	}
	else
		glDrawElements(GL_TRIANGLES, num_indices, GL_UNSIGNED_INT, nullptr);
	vao.disable(ctx);
}


GPUgeometry tessellator::quad (
	cgv::render::context &ctx, cgv::render::shader_program &shader,const vec3 &min_corner,
	const vec3 &max_corner, unsigned num_samples0, unsigned num_samples1
){
	// determine "flat" dimension of 2D-quad embedded in 3D
	const vec3 ext = max_corner - min_corner;
	const unsigned dim_flat = ext.x() < ext.y() ?   (ext.x() < ext.z() ? 0 : 2)
	                                              : (ext.y() < ext.z() ? 1 : 2);

	// determine dimension of maximal extent (will decide "major-ness" of vertex definition)
	auto [dim_i, dim_j] = ordered_unmasked_dims(ext, dim_flat);

	// create quad geometry defined by the given min/max corner and sample counts
	// - prelude
	const float i_denom = (float)num_samples0-1, j_denom = (float)num_samples1-1;
	std::vector<vertex> vertices;
	std::vector<unsigned> indices;
	GPUgeometry geom; geom.restart_index = (unsigned)-1;
	// - create actual geometry
	for (unsigned j=0; j<num_samples1; j++)
	{ 
		const float j_lerp = j/j_denom;
		for (unsigned i=0; i<num_samples0; i++)
		{
			const float i_lerp = i/i_denom;

			// create new vertex
			vertex new_vert;
			// - set attributes
			new_vert.position(dim_i) = lerp(min_corner(dim_i), max_corner(dim_i), i_lerp);
			new_vert.position(dim_j) = lerp(min_corner(dim_j), max_corner(dim_j), j_lerp);
			new_vert.position(dim_flat) = lerp(min_corner(dim_flat), max_corner(dim_flat), j_lerp); // could use i_lerp as well for slightly different result
			new_vert.position.w() = 1;
			new_vert.texcoord.set(i_lerp, j_lerp);
			// - commit
			vertices.emplace_back(new_vert);

			// determine triangle strip vertex indices for this "row"
			if (j < num_samples1-1) {
				indices.emplace_back(  j  *num_samples0 + i);
				indices.emplace_back((j+1)*num_samples0 + i);
			}
		}

		// add restart index
		if (j < num_samples1-2)
			indices.emplace_back(geom.restart_index);
	}

	// upload geometry to GPU
	// - upload attributes
	bool success = geom.vbo.create(ctx, vertices);
	success &= geom.indices.create(ctx, indices);
	// - create layout description
	success &= geom.vao.create(ctx);
	success &= geom.vao.bind_attribute_array(
		ctx, shader, "position", vertex::vec4desc, geom.vbo,
		/* offset from start of vertex buffer */ vertex::position_offset,
		/* number of elements in the array */    vertices.size(),
		/* stride between successive elements */ sizeof(vertex)
	);
	success &= geom.vao.bind_attribute_array(
		ctx, shader, "texcoord", vertex::vec2desc, geom.vbo,
		/* offset from start of vertex buffer */ vertex::texcoord_offset,
		/* number of elements in the array */    vertices.size(),
		/* stride between successive elements */ sizeof(vertex)
	);
	success &= geom.vao.set_element_array(ctx, geom.indices);

	// done!
	return success ? std::move(geom) : GPUgeometry{};
}
