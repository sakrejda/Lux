add_library(Lux SHARED survival_FLAT.cpp proposals_FLAT.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES survival_FLAT.hpp proposals_FLAT.hpp DESTINATION include)
