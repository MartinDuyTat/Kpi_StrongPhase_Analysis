// Martin Duy Tat 13th May 2021
/**
 * DoubleTagMeasurement is a class that stores the hadronic parameters and the double tag yields of a particular mode
 */

#ifndef DOUBLETAGMEASUREMENT
#define DOUBLETAGMEASUREMENT

#include<string>
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
    double GetChi2(double Normalization, double rDcosDelta, double rDsinDelta);
    /**
     * Get the number of bins
     */
    int GetNBins() const ;
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
};

#endif
