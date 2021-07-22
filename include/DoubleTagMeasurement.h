// Martin Duy Tat 13th May 2021
/**
 * DoubleTagMeasurement is a class that stores the hadronic parameters and the double tag yields of a particular mode
 */

#ifndef DOUBLETAGMEASUREMENT
#define DOUBLETAGMEASUREMENT

#include<string>
#include<vector>
#include"HadronicParameters.h"
#include"DoubleTagYields.h"

class DoubleTagMeasurement {
  public:
    /**
     * Constructor that reads in all the measurements from text files
     * @param NBins Number of bins in binning scheme
     * @param K0Mode "KS" or "KL"
     * @param HadronicParametersFilename Text file with measured hadronic parameters
     * @param DTYieldsFilename Text file with measured double tag yields
     */
    DoubleTagMeasurement(int NBins, const std::string &K0Mode, const std::string &HadronicParametersFilename, const std::string &DTYieldsFilename);
    /**
     * Function for obtaining the \f$\chi^2\f$ of the double tag yield, compared with the prediction from hadronic parameters
     */
    double GetChi2(double Normalization, double rDcosDelta, double rDsinDelta, const std::string &ErrorCategory);
    /**
     * Get the number of bins
     */
    int GetNBins() const;
    /**
     * Smear \f$K_i\f$ for systematics studies
     */
    void SmearKi();
    /**
     * Smear \f$c_i, s_i\f$ for systematics studies
     * @param A vector with smearing in the order c_1, c_2, ..., s_1, s_2, ...
     */
    void Smearcisi(const std::vector<double> &Smearing);
    /**
     * Remove all smearing
     */
    void RemoveSmearing();
    /**
     * Get the tag mode, "KSKK", "KLKK", "KSpipi" or "KLpipi"
     */
    std::string Mode() const;
  private:
    /**
     * Number of bins (not counting negative bins)
     */
    int m_NBins;
    /**
     * Hadronic parameters
     */
    HadronicParameters m_HParameters;
    /**
     * Double tag yields
     */
    DoubleTagYields m_DTYields;
    /**
     * The tag mode of this dataset
     */
    std::string m_Mode;
};

#endif
