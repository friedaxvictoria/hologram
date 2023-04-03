////
// This file should do generic, module-wide initializations (usually, those are very few)

// register the embedded shader sources with the Framework when building a single executable
#ifdef REGISTER_SHADER_FILES
	#include <cgv/base/register.h>
	#include <multi_view_rendering_shader_inc.h>
#endif
