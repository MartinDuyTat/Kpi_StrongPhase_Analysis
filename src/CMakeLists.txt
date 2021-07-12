add_library(KpiStrongPhase Analyse.cpp AnalyseBinMigration.cpp AnalyseFinalYields.cpp AnalysePeakingBackgrounds.cpp AnalyseTruthYield.cpp AnalyseYield.cpp BinningScheme.cpp Chi2DoubleTagYield.cpp DeltaKpiFitter.cpp DoubleTagMeasurement.cpp DoubleTagYields.cpp FlavourCorrection.cpp HadronicParameters.cpp KSKKBackground.cpp Settings.cpp TreeWrapper.cpp)

target_include_directories(KpiStrongPhase PUBLIC ../include)

target_link_libraries(KpiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad ROOT::Minuit2)
