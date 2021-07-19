// Martin Duy Tat 13th May 2021

#ifndef HADRONICPARAMETERS
#define HADRONICPARAMETERS

#include"BinVector.h"
#include"TMath.h"

/**
 * HadronicParameters is a class that stores the hadronic parameters \f$K_i\f$, \f$c_i\f$ and \f$s_i\f$ for a particular decay
 */
class HadronicParameters {
  public:
    /**
     * Default constructor
     */
    HadronicParameters();
    /**
     * Constructor that reads the hadronic parameters from a file
     * @param Filename Name of file with hadronic parameters
     * @param NBins Number of bins (not counting negative bins)
     * @param K0Mode "KS" for \f$K_S^0\f$ and "KL" for \f$K_L\f$
     */
    HadronicParameters(const std::string &Filename, int NBins, const std::string &K0Mode);
    /**
     * Function that creates a formula string for the predicted number of events in each bin
     * The order is: Normalization H, rD, cos(delta) and sin(delta)
     * @param Bin Bin number
     */
    std::string CreateYieldFormula(int Bin) const;
    /**
     * Function for getting \f$c_i\f$
     * @param Bin Bin number
     */
    double Getci(int Bin) const;
    /**
     * Function for getting \f$s_i\f$
     * @param Bin Bin number
     */
    double Getsi(int Bin) const;
    /**
     * Function for getting \f$K_i\f$
     * @param Bin Bin number
     */
    double GetKi(int Bin) const;
    /**
     * Function for getting the error on \f$K_i\f$
     * @param Bin Bin number
     */
    double GetKiError(int Bin) const;
    /**
     * Function that normalizes the predicted yield
     */
    double NormalizeYield(double rDcosDelta, double rDsinDelta) const;
    /**
     * Function that calculates the raw predicted yield
     */
    double CalculateYield(int Bin, double rDcosDelta, double rDsinDelta) const;
    /**
     * Function that calculates the error on the raw predicted yield from the Ki
     */
    double CalculateYieldKiError(int Bin, double rDcosDelta, double rDsinDelta) const;
    /**
     * Function that calculates the error on the raw predicted yield from the ci and si
     */
    double CalculateYieldcisiError(int Bin, double rDcosDelta, double rDsinDelta) const;
    /**
     * Function that calculates all the normalized yields and their error, and returns by reference
     * @param Yield Vector of normalized yields returned by reference
     * @param YieldError Vector of normalized yield errors returned by reference
     * @param Normalization A normalization constant for the predicted yield that is equal to 1 when the yields are normalized
     */
    void CalculateNormalizedYields(double Normalization, double rDcosDelta, double rDsinDelta, BinVector<double> &Yield, BinVector<double> &YieldKiError, BinVector<double> &YieldcisiError) const;
    /**
     * Function that prints the hadronic parameters
     */
    void PrintHadronicParameters() const;
  private:
    /**
     * Number of bins (not counting negative bins);
     */
    int m_NBins;
    /**
     * "KS" or "KL"
     */
    std::string m_K0Mode;
    /**
     * Fractional yields \f$K_i\f$
     */
    BinVector<double> m_Ki;
    /**
     * Error on the fractional yields \f$K_i\f$
     */
    BinVector<double> m_KiError;
    /**
     * Cosine of strong phase \f$c_i\f$
     */
    BinVector<double> m_ci;
    /**
     * Error on the cosine of strong phase \f$c_i\f$
     */
    BinVector<double> m_ciError;
    /**
     * Sine of strong phase \f$c_i\f$
     */
    BinVector<double> m_si;
    /**
     * Error on sine of strong phase \f$c_i\f$
     */
    BinVector<double> m_siError;
};

#endif
