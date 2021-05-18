add_library(KpiStrongPhase Chi2DoubleTagYield.cpp DeltaKpiFitter.cpp DoubleTagMeasurement.cpp DoubleTagYields.cpp HadronicParameters.cpp Settings.cpp)

target_include_directories(KpiStrongPhase PUBLIC ../include)

target_link_libraries(KpiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad ROOT::Minuit2)
