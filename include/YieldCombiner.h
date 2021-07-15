// Martin Duy Tat 14th July 2021
/**
 * YieldCombiner takes a set of binned yields, adds them, propagates the errors and normalizes the final yields
 */

#ifndef YIELDCOMBINER
#define YIELDCOMBINER

#include<string>
#include<vector>
#include"BinVector.h"

class YieldCombiner {
  public:
    /**
     * Constructor that takes a list of filenames with yields of different tags
     * @param Bins Number of bins
     * @param Files List of filenames
     */
    YieldCombiner(int Bins, const std::vector<std::string> &Files);
    /**
     * Save the final yields
     * @param Filename Very obvious!
     */
    void SaveYields(const std::string &Filename) const;
  private:
    /**
     * Binned yield
     */
    BinVector<double> m_Yield;
    /**
     * Binned yield total error, after combining statistical and systematic errors
     */
    BinVector<double> m_YieldError;
    /**
     * Function that normalizes the yields, including propagating the uncertainties
     */
    void NormalizeYields();
    /**
     * Since the normalization constant depends on the individual yields, the error propagation is more complex
     * @param Bin number to calculate error for
     * @param Sum Normalization sum of all yields
     * @param NormalizedYield The yield after normalization
     * @param Error The yield error before normalization
     */
    double CalculateNormalizationError(int Bin, double Sum, const BinVector<double> &NormalizedYield, const BinVector<double> &Error) const;
};

#endif
