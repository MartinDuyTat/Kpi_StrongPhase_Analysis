add_library(KpiStrongPhase Analyse.cpp AnalyseBinMigration.cpp AnalyseFinalYields.cpp AnalysePeakingBackgrounds.cpp AnalyseTruthYield.cpp AnalyseYield.cpp BinningScheme.cpp Chi2DoubleTagYield.cpp cisiCovariance.cpp DeltaKpiFitter.cpp DoubleTagMeasurement.cpp DoubleTagYields.cpp FlavourCorrection.cpp HadronicParameters.cpp KSKKBackground.cpp PresentYields.cpp Settings.cpp TreeWrapper.cpp YieldCombiner.cpp)

target_include_directories(KpiStrongPhase PUBLIC ../include)

target_link_libraries(KpiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad ROOT::Minuit2)
