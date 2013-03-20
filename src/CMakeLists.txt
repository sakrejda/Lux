add_library(Lux SHARED recapture.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES slicer.hpp recapture.hpp DESTINATION include)
