add_library(Lux SHARED recapture.cpp location.cpp t_walk.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES recapture.hpp location.hpp t_wlak.hpp DESTINATION include)
