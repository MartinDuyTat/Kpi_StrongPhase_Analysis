// Martin Duy Tat 19th May 2021
/**
 * AnalyseFinalYields is a class that puts together all the generator, signal MC and data yields and obtains a normalized, efficiency corrected yield in each bin
*/

#ifndef ANALYSEFINALYIELDS
#define ANALYSEFINALYIELDS

#include<string>
#include"Analyse.h"
#include"BinVector.h"

class AnalyseFinalYields: public Analyse {
  public:
    /**
     * Constructor that loads the generator, signal MC and data yields
     * @param GeneratorYieldsFilename Filename of text file with generator yields
     * @param SignalMCYieldsFilename Filename of text file with signal MC yields
     * @param DataYieldsFilename Filename of text file with data yields
     * @param FlavourTagCorrectionFilename Filename of text file with flavour tag corrections
     */
    AnalyseFinalYields(const std::string &GeneratorYieldsFilename, const std::string &SignalMCYieldsFilename, const std::string &DataYieldsFilename, const std::string &FlavourTagCorrectionFilename = "");
    /**
     * Function that calculates the efficiency corrected, normalized yields and saves them to a file
     * @param Filename Filename of text file with final yields
     */
    void CalculateFinalYields(const std::string &Filename) const;
  private:
    /**
     * Vector with generator level yields
     */
    BinVector<double> m_GeneratorYields;
    /**
     * Vector with signal MC yields
     */
    BinVector<double> m_SignalMCYields;
    /**
     * Vector with data yields
     */
    BinVector<double> m_DataYields;
    /**
     * Vector with data yield errors
     */
    BinVector<double> m_DataYieldErrors;
    /**
     * Flavour tag corrections
     */
    BinVector<double> m_FlavourTagCorrections;
    /**
     * Flavour tag correction errors
     */
    BinVector<double> m_FlavourTagCorrectionErrors;
};

#endif
