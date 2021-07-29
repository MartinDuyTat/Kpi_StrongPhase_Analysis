// Martin Duy Tat 13th May 2021

#ifndef CHI2DOUBLETAGYIELD
#define CHI2DOUBLETAGYIELD

#include<vector>
#include<string>
#include<utility>
#include"DoubleTagMeasurement.h"
#include"cisiCovariance.h"
#include"Minuit2/Minuit2Minimizer.h"

/**
 * Chi2DoubleTagYield is a class for combining different measurements and getting the total \f$\chi^2\f$
*/
class Chi2DoubleTagYield {
  public:
    /**
     * Constructor that fixed the normalization by default
     * @param FixNormalization True if normalization is fixed to 1
     * @param ErrorCategory The type of error considered in the fit, "Kpi", "Ki" or "cisi"
     */
    Chi2DoubleTagYield(bool FixNormalization, const std::string &ErrorCategory = "Kpi");
    /**
     * Function for adding a measurement
     * @param NBins Number of bins
     * @param K0Mode "KSKK", "KLKK", "KSpipi", or "KLpipi"
     * @param DataSetsToFit String containing (any) K0Mode that should be included in the fit
     * @param cisiHParameterFilename File with strong phase \f$c_i, s_i\f$ parameters
     * @param KiHParameterFilename File with hadronic \f$K_i\f$ parameters
     * @param DTYieldFilename File with double tag yields
     */
    void AddMeasurement(int NBins, const std::string &K0Mode, const std::string &DataSetsToFit, const std::string &cisiHParameterFilename, const std::string &KiHParameterFilename, const std::string &DTYieldFilename);
    /**
     * () operator overload to easily get the total \f$\chi^2\f$
     */
    double operator()(const double *params);
    /**
     * Run Minuit to minimize the \f$\chi^2\f$
     */
    void MinimizeChi2(const std::string &PlotContourFilename = "None");
    /**
     * Function that plots the 1, 2, 3 sigma contours of the results
     */
    void DrawContours(ROOT::Minuit2::Minuit2Minimizer *Minimizer, const std::string &Filename) const;
    /**
     * Get fitted value of \f$r_D\cos(\delta)\f$
     */
    double GetFittedrDcosDelta() const;
    /**
     * Get fitted error of \f$r_D\cos(\delta)\f$
     */
    double GetErrorrDcosDelta() const;
    /**
     * Get fitted value of \f$r_D\sin(\delta)\f$
     */
    double GetFittedrDsinDelta() const;
    /**
     * Get fitted value of \f$r_D\sin(\delta)\f$
     */
    double GetErrorrDsinDelta() const;
    /**
     * Get the correlation coefficient between \f$r_D\cos(\delta)\f$ and \f$r_D\sin(\delta)\f$
     */
    double GetCorrelation() const;
    /**
     * Get the $\chi^2\f$ per degrees of freedom
     */
    double GetChi2PerDegreesOfFreedom() const;
    /**
     * Get the number of degree of freedom
     */
    int GetDegreesOfFreedom() const;
    /**
     * Run \f$K_i\f$ systematics by smearing \f$K_i\f$ and do 1000 fits
     * @param Systematics "Ki" or "cisi"
     * @param rDcosDelta_Bias The bias in \f$r_D\cos(\delta)\f$ under smearing
     * @param rDsinDelta_Bias The bias in \f$r_D\sin(\delta)\f$ under smearing
     * @param rDcosDelta_Syst The systematic uncertainty in \f$r_D\cos(\delta)\f$ under smearing
     * @param rDsinDelta_Syst The systematic uncertainty in \f$r_D\sin(\delta)\f$ under smearing
     * @param Correlation The systematic correlation between \f$r_D\cos(\delta)\f$ and \f$r_D\sin(\delta)\f$
     */
    void RunSystematics(const std::string &Systematics, double &rDcosDelta_Bias, double &rDsinDelta_Bias, double &rDcosDelta_Syst, double &rDsinDelta_Syst, double &Correlation);
    /**
     * Setter for m_VetoBins
     */
    void SetVetoBins(const std::vector<std::pair<std::string, int>> &VetoBins);
  private:
    /**
     * Flag that is true when the normalization is fixed to 1
     */
    bool m_FixNormalization;
    /**
     * The error type considered, "Kpi", "Ki" or "cisi"
     */
    std::string m_ErrorCategory;
    /**
     * Vector of measurements
     */
    std::vector<DoubleTagMeasurement> m_Measurements;
    /**
     * Fitted value of \f$r_D\cos(\delta)\f$
     */
    double m_FittedrDcosDelta;
    /**
     * Fitted error of \f$r_D\cos(\delta)\f$
     */
    double m_ErrorrDcosDelta;
    /**
     * Fitted value of \f$r_D\sin(\delta)\f$
     */
    double m_FittedrDsinDelta;
    /**
     * Fitted value of \f$r_D\sin(\delta)\f$
     */
    double m_ErrorrDsinDelta;
    /**
     * The correlation coefficient between \f$r_D\cos(\delta)\f$ and \f$r_D\sin(\delta)\f$
     */
    double m_Correlation;
    /**
     * \f$\chi^2\f$ per degree of freedom of fit
     */
    double m_Chi2;
    /**
     * Object containing all the covariance matrices for the \f$c_i, s_i\f$
     */
    cisiCovariance m_cisiCovariance;
    /**
     * Vector of pairs, the first component being the K0 mode and the other being the bin that is left out from the fit
     */
    std::vector<std::pair<std::string, int>> m_VetoBins;
    /**
     * Helper function that sets up Minuit2 and runs the minimization
     */
    ROOT::Minuit2::Minuit2Minimizer* RunMinimization(double &rDcosDelta, double &rDsinDelta, double &rDcosDeltaError, double &rDsinDeltaError, double &Chi2, double &Correlation) const;
};

#endif
