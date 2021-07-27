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
     * @param DataSetsToFit A string containing all the K0hh modes that should be added in the fit
     * @param FixNormalization Flag that fixes the normalization if it's true
     * @param ErrorCategory The type of error considered in the fit, "Kpi", "Ki" or "cisi"
     * @param VetoBinsFilename Filename of text file containing modes and bins that are left out of the fit
     */
    DeltaKpiFitter(const std::string &Filename, const std::string &DataSetsToFit, bool FixNormalization = true, const std::string &ErrorCategory = "Kpi", const std::string &VetoBinsFilename = "None");
    /**
     * Run the fit, run \f$K_i\f$ systematics, run \f$c_i, s_i\f$ systematics and save measurements
     * @param Filename Filename of text file to save results to
     * @param PlotContourFilename Filename for contours, leave empty to not plot them
     * @param RunKiSystematics Set to true to do a \f$K_i\f$ systematics study
     * @param RuncisiSystematics Set to true to do a \f$c_i, s_i\f$  systematics study
     */
    void RunFit(const std::string &Filename, const std::string &PlotContourFilename = "None", bool RunKiSystematics = false, bool RuncisiSystematics = false);
  private:
    /**
     * Object containing the measurements
     */
    Chi2DoubleTagYield m_Measurements;
};

#endif
