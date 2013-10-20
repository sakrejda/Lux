add_library(Lux SHARED 
	recapture.cpp 
	locations.cpp 
	rv_constant.cpp
	rv_uniform.cpp
	rv_normal.cpp
	rv_t_walk.cpp
	rv_missing_t_walk.cpp
	rv_missing_t_walk_observed_normal.cpp
	rv_missing_t_walk_observed_interval.cpp
)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES recapture.hpp locations.hpp t_walk.hpp random.hpp DESTINATION include)
