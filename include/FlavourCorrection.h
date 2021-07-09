// Martin Duy Tat 9th July 2021

#ifndef FLAVOURCORRECTION
#define FLAVOURCORRECTION

#include<string>
#include<map>
#include<utility>
#include"BinVector.h"

/**
 * FlavourCorrection is a class that calculates the DCS flavour tag corrections based on amplitude-averaged strong phases and hadronic parameters
 */

class FlavourCorrection {
  public:
    /**
     * Constructor that takes in files with ci, si, Ki and hadronic parameters
     * @param Bins Number of bins
     * @param Mode "KSKK" or "KLKK"
     * @param StrongPhaseFilename Filename of file with ci, si, Ki
     * @param HadronicFilename Filename with rD, deltaD and R
     */
    FlavourCorrection(int Bins, const std::string &Mode, const std::string &StrongPhaseFilename, const std::string &HadronicFilename);
    /**
     * Calculate the flavour tag correction in bin i and its error
     */
    std::pair<double, double> CalculateCorrection(int Bin) const;
    /**
     * Save corrections to file
     * @param Filename File to save corrections to
     */
    void SaveCorrections(const std::string &Filename) const;
  private:
    /**
     * Number of bins
     */
    int m_Bins;
    /**
     * The mode "KSKK" or "KLKK"
     */
    std::string m_Mode;
    /**
     * The amplitude-averaged cosine of the strong phase
     */
    BinVector<double> m_ci;
    /**
     * The amplitude-averaged sine of the strong phase
     */
    BinVector<double> m_si;
    /**
     * The fractional yield
     */
    BinVector<double> m_Ki;
    /**
     * The hadronic rD, deltaD and R parameters
     */
    std::map<std::string, double> m_Hadronic;
    /**
     * The hadronic rD, deltaD and R parameter errors
     */
    std::map<std::string, double> m_HadronicErrors;
};

#endif
