// Martin Duy Tat 22nd July 2021
/**
 * cisiCovariance is a class that stores the covariance matrices of KSKK, KLKK, KSpipi, KLpipi tag mode \f$c_i, s_i\f$, both statistical and systematic, calculates the Cholesky decomposition of these and then uses them to generate correlated smearing of the strong phases
 */

#ifndef CISICOVARIANCE
#define CISICOVARIANCE

#include<string>
#include<vector>
#include<map>
#include"BinVector.h"
#include"TMatrixT.h"

class cisiCovariance {
  public:
    /**
     * Constructor that initializes all the vectors and matrices to the correct size
     */
    cisiCovariance();
    /**
     * Read the uncertainties and correlation matrices from a file
     * It is assumed that the correlation matrix is read from the file with the KShh mode, only uncertainties are read from the KLhh file
     */
    void AddDataset(int NBins, const std::string &Mode, const std::string &cisiHadronicParameterFilename);
    /**
     * Once all the correlation matrices have been read from file, this function will perform the Cholesky decomposition
     * If not all the tag modes have been read from file yet, an exception is thrown
     */
    void PrepareCholesky();
    /**
     * Generate the smearing for all \f$c_i, s_i\f$ for all datasets
     */
    void Smear();
    /**
     * Get the correct smearing
     * @param Mode "KSKK", "KLKK", "KSpipi" or "KLpipi"
     */
    std::vector<double> GetSmearing(const std::string &Mode) const;
  private:
    /**
     * A map of all tag modes (for looping) and a flag indicating if their covariance matrix has been filled
     */
    std::map<std::string, bool> m_Modes;
    /**
     * A map of statistical correlation matrices for the different tag modes
     */
    std::map<std::string, TMatrixT<double>> m_StatCorr;
    /**
     * A map of systematic correlation matrices for the different tag modes
     */
    std::map<std::string, TMatrixT<double>> m_SystCorr;
    /**
     * A map of the statistical uncertainties for all the \f$c_i, s_i\f$ for the different tag modes
     */
    std::map<std::string, std::vector<double>> m_StatUncertainties;
    /**
     * A map of the systematic uncertainties for all the \f$c_i, s_i\f$ for the different tag modes
     */
    std::map<std::string, std::vector<double>> m_SystUncertainties;
    /**
     * A map of the statistical Cholesky matrix for the different tag modes
     */
    std::map<std::string, TMatrixT<double>> m_StatCholesky;
    /**
     * A map of the systematic Cholesky matrix for the different tag modes
     */
    std::map<std::string, TMatrixT<double>> m_SystCholesky;
    /**
     * A map of the smearing of all \f$c_i, s_i\f$ for the different tag modes
     */
    std::map<std::string, std::vector<double>> m_Smearing;
};

#endif
