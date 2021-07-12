// Martin Duy Tat 12th July 2021
/**
 * KSKKBackground is a class for calculating the binned yields of KSKK background in KLKK tags
 * The prediction is based on the total background yield, the strong phases of KSKK from data, fractional yields from the model and the hadronic parameters \f$r_D\f$ and \f$\delta_D\f$
 */

#ifndef KSKKBACKGROUND
#define KSKKBACKGROUND

#include<vector>
#include"BinVector.h"

class KSKKBackground {
  public:
    /**
     * Constructor that sets up all the binned yields and their errors
     */
    KSKKBackground();
    /**
     * Function that returns the binned yield
     * @param Bin
     */
    double GetBinYield(int Bin) const;
    /**
     * Function that returns the binned yield errorn
     * @param Bin
     */
    double GetBinYieldError(int Bin) const;
  private:
    /**
     * Vector with binned yields
     */
    BinVector<double> m_Yields;
    /**
     * Vector with binned yield errors
     */
    BinVector<double> m_YieldErrors;
};

#endif
