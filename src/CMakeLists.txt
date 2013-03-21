add_library(Lux SHARED recapture.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES recapture.hpp DESTINATION include)
