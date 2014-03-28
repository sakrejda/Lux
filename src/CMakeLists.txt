add_library(Lux SHARED 
	recapture.cpp 
	time_series.cpp 
	rv_constant.cpp
	rv_uniform.cpp
	rv_normal.cpp
	rv_t_walk.cpp
	rv_missing_t_walk_core.cpp
	rv_missing_t_walk_observed_normal.cpp
	rv_missing_t_walk_observed_interval.cpp
	rv_missing_t_walk.cpp
)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES 
	recapture.hpp 
	time_series.hpp
	random.hpp 
	rv_constant.hpp
	rv_uniform.hpp
	rv_normal.hpp
	rv_t_walk.hpp
	rv_missing_t_walk_core.hpp
	rv_missing_t_walk_observed_normal.hpp
	rv_missing_t_walk_observed_interval.hpp
	rv_missing_t_walk.hpp
DESTINATION include)
