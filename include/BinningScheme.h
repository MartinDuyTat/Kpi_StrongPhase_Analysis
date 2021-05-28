// Martin Duy Tat 18th May 2021
/**
 * BinningScheme is class that contains the binning scheme, and one can determine the bin number from the Dalitz coordinates
 */

#ifndef BINNINGSCHEME
#define BINNINGSCHEME

#include"TH2F.h"

class BinningScheme {
  public:
    /**
     * Constructor that loads the binning scheme into memory
     */
    BinningScheme();
    /**
     * Get the bin number
     * @param M2Plus Invariant mass of \f$K_S^0\f$ and \f$K^+\f$
     * @param M2Plus Invariant mass of \f$K_S^0\f$ and \f$K^-\f$
     * @param KCharge Charge of the tag kaon
     */
    int GetBinNumber(double M2Plus, double M2Minus, int KCharge) const;
    /**
     * Get the closest bin number, assuming the initial point is outside the Dalitz plot
     * Algorithm by Sophia Vaughan
     * @param M2Plus Invariant mass of \f$K_S^0\f$ and \f$K^+\f$
     * @param M2Plus Invariant mass of \f$K_S^0\f$ and \f$K^-\f$
     * @param KCharge Charge of the tag kaon
     */
    int GetMappedBinNumber(double M2Plus, double M2Minus, int KCharge) const;
    /**
     * Get number of bins
     */
    int GetNumberBins() const;
  private:
    /**
     * 2D histogram that contains the binning scheme
     */
    TH2F *m_BinningScheme;
    /**
     * Number of bins in binning scheme
     */
    int m_NumberBins;
};

#endif
