#pragma once

#include <cgv/math/fvec.h>
#include <cgv/math/fmat.h>
#include <cgv/media/axis_aligned_box.h>
#include "lib_begin.h"


////
// misc utilities

/// de-homogenizes the given vector by dividing each component by the last component and then removing the latter
template <class T, unsigned N>
cgv::math::fvec<T, N-1> w_clip (const cgv::math::fvec<T, N> &v)
{
	const T inv_w = 1/v(N-1);
	cgv::math::fvec<T, N-1> result;
	for (unsigned i=0; i<N-1; i++)
		result(i) = v(i) * inv_w;
	return result;
}

/// finds the near clipping plane from the given projection matrix
/// NOTE: not optimized - calculates the inverse of the input matrix, could be done in different more efficient ways
template <class T>
T znear_from_mat (const cgv::math::fmat<T, 4, 4> &projection_matrix)
{
	static cgv::math::fvec<T, 3> p_near(0, 0, 0);
	const auto inv_proj = cgv::math::inv(projection_matrix);
	const auto zn_vec = inv_proj.mul_pos(p_near);
	return -zn_vec.z() / zn_vec.w();
}

/// finds the far clipping plane from the given projection matrix
/// NOTE: not optimized - calculates the inverse of the input matrix, could be done in different more efficient ways
template <class T>
T zfar_from_mat (const cgv::math::fmat<T, 4, 4> &projection_matrix)
{
	static cgv::math::fvec<T, 3> p_far(0, 0, 1);
	const auto inv_proj = cgv::math::inv(projection_matrix);
	const auto zf_vec = inv_proj.mul_pos(p_far);
	return -zf_vec.z() / zf_vec.w();
}

/// finds the projection of the given box onto the given direction vector
template <class T>
const std::pair<const T, const T> project_box_onto_dir (
	const cgv::media::axis_aligned_box<T, 3> &box, const cgv::math::fvec<T, 3> d
)
{
	const auto he = box.get_extent()/2;
	const auto center = box.get_center();
	const T center_p = cgv::math::dot(d, center);
	const T radius = he[0]*std::abs(d[0]) + he[1]*std::abs(d[1]) + he[2]*std::abs(d[2]);
	return {center_p-radius, center_p+radius};
}
