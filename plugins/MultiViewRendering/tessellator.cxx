#include <memory>
#include <libs/cgv_gl/gl/gl.h>
#include <cgv/math/functions.h>
#include "tessellator.h"

using namespace cgv::math;
using namespace cgv::render;

// anonymous namespace for helper functions
namespace {
// given a "mask" number m∈{0,1,2}, returns the other two (i.e. unmasked) numbers
inline std::pair<unsigned, unsigned> select_unmasked(const unsigned mask)
{
	return {mask != 0 ? 0 : (mask != 1 ? 1 : 2), mask != 2 ? 2 : 1};
}

// given a number i∈{0,1,2} "masking" a dimension from the given 3-vector, return the indices
// of the other dimensions ordered according to the value of their corresponding components
// (index of larger component first, index of smaller component second)
template <class T>
inline const std::pair<const unsigned, const unsigned> ordered_unmasked_dims(const cgv::math::fvec<T, 3>& vec,
																			 const unsigned mask)
{
	const auto& [dim0, dim1] = select_unmasked(mask);
	if (vec(dim1) > vec(dim0))
		return {dim1, dim0};
	else
		return {dim0, dim1};
}

// helper to determine the number of required float components for a given attribute combination
constexpr unsigned num_components(const tessellator::vertex_attribs attribs)
{
	switch (attribs) {
	case tessellator::VA_POSITION:
		return 4;
	case tessellator::VA_NORMAL:
		return 3;
	case tessellator::VA_POS_NORMAL:
		return num_components(tessellator::VA_POSITION) + num_components(tessellator::VA_NORMAL);
	case tessellator::VA_TEXCOORD:
		return 2;
	case tessellator::VA_POS_TEXCOORD:
		return num_components(tessellator::VA_POSITION) + num_components(tessellator::VA_TEXCOORD);
	case tessellator::VA_NORMAL_TEXCOORD:
		return num_components(tessellator::VA_NORMAL) + num_components(tessellator::VA_TEXCOORD);
	case tessellator::VA_POS_NORMAL_TEXCOORD:
		return num_components(tessellator::VA_POSITION) + num_components(tessellator::VA_NORMAL) +
			   num_components(tessellator::VA_TEXCOORD);
	default:
		return -1;
	}
}

// the format of our vertices
struct vertex_store
{
	typedef typename tessellator::vec2 vec2;
	typedef typename tessellator::vec3 vec3;
	typedef typename tessellator::vec4 vec4;

	vertex_store(const tessellator::vertex_attribs attribs, const unsigned num)
		: attribs(attribs), num(num), num_floats(num_components(attribs) * num), size(num_floats * sizeof(float)),
		  stride(num_components(attribs) * sizeof(float)), c(std::make_unique<float[]>(num_floats))
	{
	}

	// the vertex attributes in this store
	const tessellator::vertex_attribs attribs;

	// the number of vertices in the store
	const unsigned num;

	// the number of float components of all vertices in the store
	const unsigned num_floats;

	// the size in bytes of all vertice in the store
	const size_t size;

	// the stride in bytes from the start of one vertex to the next
	const unsigned stride;

	// the actual vertex component store
	std::unique_ptr<float[]> c;

	inline static constexpr unsigned pos_offset(const tessellator::vertex_attribs attribs)
	{
		return attribs & tessellator::VA_POSITION ? 0 : -1;
	}
	inline static constexpr unsigned normal_offset(const tessellator::vertex_attribs attribs)
	{
		return attribs & tessellator::VA_NORMAL
					 ? (attribs & tessellator::VA_POSITION ? num_components(tessellator::VA_POSITION) : 0)
					 : -1;
	}
	inline static constexpr unsigned texcoord_offset(const tessellator::vertex_attribs attribs)
	{
		return attribs & tessellator::VA_TEXCOORD
					 ? (attribs & tessellator::VA_POSITION
							  ? (attribs & tessellator::VA_NORMAL ? num_components(tessellator::VA_POS_NORMAL)
																  : num_components(tessellator::VA_POSITION))
							  : 0)
					 : -1;
	}
	inline unsigned pos_offset(void) const { return pos_offset(attribs); }
	inline unsigned normal_offset(void) const { return normal_offset(attribs); }
	inline unsigned texcoord_offset(void) const { return normal_offset(attribs); }

	inline static constexpr unsigned pos_offset_bytes(const tessellator::vertex_attribs attribs)
	{
		return pos_offset(attribs) * sizeof(float);
	}
	inline static constexpr unsigned normal_offset_bytes(const tessellator::vertex_attribs attribs)
	{
		return normal_offset(attribs) * sizeof(float);
	}
	inline static constexpr unsigned texcoord_offset_bytes(const tessellator::vertex_attribs attribs)
	{
		return texcoord_offset(attribs) * sizeof(float);
	}
	inline unsigned pos_offset_bytes(void) const { return pos_offset_bytes(attribs); }
	inline unsigned normal_offset_bytes(void) const { return normal_offset_bytes(attribs); }
	inline unsigned texcoord_offset_bytes(void) const { return texcoord_offset_bytes(attribs); }

	inline vec4& pos(const unsigned vertex) { return *(vec4*)(c.get() + num_components(attribs) * vertex); }
	inline vec3& normal(const unsigned vertex)
	{
		return *(vec3*)(c.get() + num_components(attribs) * vertex + normal_offset(attribs));
	}
	inline vec2& texcoord(const unsigned vertex)
	{
		return *(vec2*)(c.get() + num_components(attribs) * vertex + texcoord_offset(attribs));
	}

	inline static cgv::render::type_descriptor vec2desc(void)
	{
		return element_descriptor_traits<vec2>::get_type_descriptor(vec2());
	}
	inline static cgv::render::type_descriptor vec3desc(void)
	{
		return element_descriptor_traits<vec3>::get_type_descriptor(vec3());
	}
	inline static cgv::render::type_descriptor vec4desc(void)
	{
		return element_descriptor_traits<vec4>::get_type_descriptor(vec4());
	}
};
}; // namespace

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

GPUgeometry tessellator::quad(cgv::render::context& ctx, cgv::render::shader_program& shader, const vec3& min_corner,
							  const vec3& max_corner, unsigned num_samples0, unsigned num_samples1,
							  vertex_attribs attribs)
{
	// determine "flat" dimension of 2D-quad embedded in 3D
	const vec3 ext = max_corner - min_corner;
	const unsigned dim_flat = ext.x() < ext.y() ? (ext.x() < ext.z() ? 0 : 2) : (ext.y() < ext.z() ? 1 : 2);

	// determine dimension of maximal extent (will decide "major-ness" of vertex definition)
	const auto& [dim_i, dim_j] = ordered_unmasked_dims(ext, dim_flat);

	// create quad geometry defined by the given min/max corner and sample counts
	// - prelude
	const float i_denom = (float)num_samples0 - 1, j_denom = (float)num_samples1 - 1;
	vertex_store verts(attribs, num_samples0 * num_samples1);
	std::vector<unsigned> indices;
	GPUgeometry geom;
	geom.restart_index = (unsigned)-1;
	// - precalculate common normal (TODO: untested!!!)
	const vec3 common_normal = [&, &dim_i = dim_i, &dim_j = dim_j]() {
		vec3 b0, b1;
		b0(dim_i) = max_corner(dim_i);
		b0(dim_j) = min_corner(dim_j);
		b0(dim_flat) = min_corner(dim_flat);
		b1(dim_i) = min_corner(dim_i);
		b1(dim_j) = max_corner(dim_j);
		b1(dim_flat) = max_corner(dim_flat);
		return normalize(cross(b0 - min_corner, b1 - min_corner));
	}();
	// - create actual geometry
	for (unsigned j = 0; j < num_samples1; j++) {
		const float j_lerp = j / j_denom;
		for (unsigned i = 0; i < num_samples0; i++) {
			const unsigned vid = j * num_samples0 + i;
			const float i_lerp = i / i_denom;

			// set attributes
			if (attribs & tessellator::VA_POSITION) {
				auto& pos = verts.pos(vid);
				pos(dim_i) = lerp(min_corner(dim_i), max_corner(dim_i), i_lerp);
				pos(dim_j) = lerp(min_corner(dim_j), max_corner(dim_j), j_lerp);
				pos(dim_flat) = lerp(min_corner(dim_flat), max_corner(dim_flat),
									 j_lerp); // could use i_lerp as well for slightly different result
				pos.w() = 1;
			}
			if (attribs & tessellator::VA_NORMAL)
				verts.normal(vid) = common_normal;
			if (attribs & tessellator::VA_TEXCOORD)
				verts.texcoord(vid).set(i_lerp, j_lerp);

			// determine triangle strip vertex indices for this "row"
			if (j < num_samples1 - 1) {
				indices.emplace_back(j * num_samples0 + i);
				indices.emplace_back((j + 1) * num_samples0 + i);
			}
		}

		// add restart index
		if (j < num_samples1 - 2)
			indices.emplace_back(geom.restart_index);
	}

	// upload geometry to GPU
	// - upload attributes
	bool success = geom.vbo.create(ctx, verts.c.get(), verts.num_floats);
	success &= geom.indices.create(ctx, indices);
	// - create layout description
	success &= geom.vao.create(ctx);
	if (attribs & tessellator::VA_POSITION)
		success &= geom.vao.bind_attribute_array(ctx, shader, "position", vertex_store::vec4desc(), geom.vbo,
												 /* offset from start of vertex buffer */ verts.pos_offset_bytes(),
												 /* number of elements in the array */ verts.num,
												 /* stride between successive elements */ verts.stride);
	if (attribs & tessellator::VA_NORMAL)
		success &= geom.vao.bind_attribute_array(ctx, shader, "normal", vertex_store::vec3desc(), geom.vbo,
												 /* offset from start of vertex buffer */ verts.normal_offset_bytes(),
												 /* number of elements in the array */ verts.num,
												 /* stride between successive elements */ verts.stride);
	if (attribs & tessellator::VA_TEXCOORD)
		success &= geom.vao.bind_attribute_array(ctx, shader, "texcoord", vertex_store::vec2desc(), geom.vbo,
												 /* offset from start of vertex buffer */ verts.texcoord_offset_bytes(),
												 /* number of elements in the array */ verts.num,
												 /* stride between successive elements */ verts.stride);
	success &= geom.vao.set_element_array(ctx, geom.indices);

	// done!
	return success ? std::move(geom) : GPUgeometry{};
}