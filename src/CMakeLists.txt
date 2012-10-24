add_library(Lux SHARED survival.cpp)
INSTALL(TARGETS Lux LIBRARY DESTINATION lib)
INSTALL(FILES survival.hpp DESTINATION include)
