#pragma once

#include <cgv_gl/renderer.h>

#include "lib_begin.h"

using namespace cgv::render;

namespace cgv {
namespace glutil {

class render_data_base : public render_types {
protected:
	attribute_array_manager aam;
	bool out_of_date;

	std::vector<unsigned> idx;
	std::vector<vec3> pos;

	bool transfer(context& ctx, renderer& r) {
		if(pos.size() > 0) {
			r.set_position_array(ctx, pos);
			if(idx.size() > 0)
				r.set_indices(ctx, idx);
			else
				r.remove_indices(ctx);
			out_of_date = false;
			return true;
		}
		return false;
	}

public:
	render_data_base() {
		out_of_date = true;
	}

	~render_data_base() { clear(); }

	bool init(context& ctx) {
		return aam.init(ctx);
	}

	void destruct(context& ctx) {
		clear();
		aam.destruct(ctx);
	}

	void clear() {
		idx.clear();
		pos.clear();
		
		out_of_date = true;
	}

	void set_out_of_date() { out_of_date = true; }

	size_t size() { return pos.size(); }

	const attribute_array_manager& ref_aam() const { return aam; }

	std::vector<unsigned>&	ref_idx() { return idx; }
	std::vector<vec3>&		ref_pos() { return pos; }

	void add_idx(const unsigned int i) {
		idx.push_back(i);
	}
};

}
}

#include <cgv/config/lib_end.h>
