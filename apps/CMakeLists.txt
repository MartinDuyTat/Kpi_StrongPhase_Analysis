add_executable(AnalyseBackgrounds AnalyseBackgrounds.cpp)
add_executable(AnalyseDoubleTagYields AnalyseDoubleTagYields.cpp)
add_executable(FitDelta FitDelta.cpp)

target_link_libraries(AnalyseBackgrounds PUBLIC KpiStrongPhase)
target_link_libraries(AnalyseBackgrounds PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(AnalyseDoubleTagYields PUBLIC KpiStrongPhase)
target_link_libraries(AnalyseDoubleTagYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitDelta PUBLIC KpiStrongPhase)
target_link_libraries(FitDelta PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS AnalyseBackgrounds AnalyseDoubleTagYields FitDelta DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
