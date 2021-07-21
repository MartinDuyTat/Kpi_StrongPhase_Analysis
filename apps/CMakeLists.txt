add_executable(AnalyseBackgrounds AnalyseBackgrounds.cpp)
add_executable(AnalyseDoubleTagYields AnalyseDoubleTagYields.cpp)
add_executable(CombineYields CombineYields.cpp)
add_executable(FitDelta FitDelta.cpp)
add_executable(DetermineFlavourTagCorrections DetermineFlavourTagCorrections.cpp)
add_executable(PlotDataWithInclusiveMC PlotDataWithInclusiveMC.cpp)
add_executable(PlotKpiYields PlotKpiYields.cpp)

target_link_libraries(AnalyseBackgrounds PUBLIC KpiStrongPhase)
target_link_libraries(AnalyseBackgrounds PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(AnalyseDoubleTagYields PUBLIC KpiStrongPhase)
target_link_libraries(AnalyseDoubleTagYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(CombineYields PUBLIC KpiStrongPhase)
target_link_libraries(CombineYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitDelta PUBLIC KpiStrongPhase)
target_link_libraries(FitDelta PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(DetermineFlavourTagCorrections PUBLIC KpiStrongPhase)
target_link_libraries(DetermineFlavourTagCorrections PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(PlotDataWithInclusiveMC PUBLIC KpiStrongPhase)
target_link_libraries(PlotDataWithInclusiveMC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(PlotKpiYields PUBLIC KpiStrongPhase)
target_link_libraries(PlotKpiYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS AnalyseBackgrounds AnalyseDoubleTagYields CombineYields FitDelta DetermineFlavourTagCorrections PlotDataWithInclusiveMC PlotKpiYields DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
