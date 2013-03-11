add_library(Lux SHARED survival_FLAT.cpp proposals_FLAT.cpp recapture.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES slicer.hpp survival_FLAT.hpp proposals_FLAT.hpp recapture.hpp DESTINATION include)
