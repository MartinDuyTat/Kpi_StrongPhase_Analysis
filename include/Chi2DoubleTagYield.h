// Martin Duy Tat 13th May 2021

#ifndef CHI2DOUBLETAGYIELD
#define CHI2DOUBLETAGYIELD

#include<vector>
#include"DoubleTagMeasurement.h"

/**
 * Chi2DoubleTagYield is a class for combining different measurements and getting the total \f$\chi^2\f$
*/
class Chi2DoubleTagYield {
  public:
    /**
     * Trivial constructor
     */
    Chi2DoubleTagYield();
    /**
     * Function for adding a measurement
     * @param Measurement A DoubleTagMeasurement object to be added
     */
    void AddMeasurement(const DoubleTagMeasurement &Measurement);
    /**
     * () operator overload to easily get the total \f$\chi^2\f$
     */
    double operator()(const double *params);
    /**
     * Run Minuit to minimiza the \f$\chi^2\f$
     */
    void MinimizeChi2();
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
     * Get the $\chi^2\f$ per degrees of freedom
     */
    double GetChi2PerDegreesOfFreedom() const;
    /**
     * Get the number of degree of freedom
     */
    int GetDegreesOfFreedom() const;
  private:
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
     * \f$\chi^2\f$ per degree of freedom of fit
     */
    double m_Chi2;
};

#endif
