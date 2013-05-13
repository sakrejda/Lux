add_library(Lux SHARED recapture.cpp location.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES recapture.hpp location.hpp DESTINATION include)
