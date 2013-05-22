add_library(Lux SHARED recapture.cpp locations.cpp t_walk.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES recapture.hpp locations.hpp t_walk.hpp DESTINATION include)
