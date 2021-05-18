add_executable(FitDelta FitDelta.cpp)

target_link_libraries(FitDelta PUBLIC KpiStrongPhase)
target_link_libraries(FitDelta PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS FitDelta DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
