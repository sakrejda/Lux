add_library(Lux SHARED recapture.cpp locations.cpp t_walk.cpp random.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES recapture.hpp locations.hpp t_walk.hpp random.hpp DESTINATION include)
