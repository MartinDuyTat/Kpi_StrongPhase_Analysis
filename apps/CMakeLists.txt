add_executable(AnalyseBackgrounds AnalyseBackgrounds.cpp)
add_executable(AnalyseDoubleTagYields AnalyseDoubleTagYields.cpp)
add_executable(FitDelta FitDelta.cpp)
add_executable(DetermineFlavourTagCorrections DetermineFlavourTagCorrections.cpp)
add_executable(PlotDataWithInclusiveMC PlotDataWithInclusiveMC.cpp)

target_link_libraries(AnalyseBackgrounds PUBLIC KpiStrongPhase)
target_link_libraries(AnalyseBackgrounds PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(AnalyseDoubleTagYields PUBLIC KpiStrongPhase)
target_link_libraries(AnalyseDoubleTagYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitDelta PUBLIC KpiStrongPhase)
target_link_libraries(FitDelta PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(DetermineFlavourTagCorrections PUBLIC KpiStrongPhase)
target_link_libraries(DetermineFlavourTagCorrections PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(PlotDataWithInclusiveMC PUBLIC KpiStrongPhase)
target_link_libraries(PlotDataWithInclusiveMC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS AnalyseBackgrounds AnalyseDoubleTagYields FitDelta DetermineFlavourTagCorrections PlotDataWithInclusiveMC DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
