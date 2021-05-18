// Martin Duy Tat 18th May 2021
/**
 * DeltaKpiFitter is a class for obtaining \f$r_D\cos(\delta)\f$ and \f$r_D\sin(\delta)\f$ from a set of measurements
 */

#ifndef DELTAKPIFITTER
#define DELTAKPIFITTER

#include<string>
#include"Chi2DoubleTagYield.h"

class DeltaKpiFitter {
  public:
    /**
     * Constructor that sets up the measurements that are fitted
     * @param Filename Filename of text file with paths to different measurements
     */
    DeltaKpiFitter(const std::string &Filename);
    /**
     * Run the fit and save measurements
     * @param Filename Filename of text file to save results to
     */
    void RunFit(const std::string &Filename);
  private:
    /**
     * Object containing the measurements
     */
    Chi2DoubleTagYield m_Measurements;
};

#endif
