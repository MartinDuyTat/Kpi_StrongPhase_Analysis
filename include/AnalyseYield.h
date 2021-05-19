// Martin Duy Tat 19th May 2021
/**
 * AnalyseYield is a class for calculating the double tag yield from reconstructed events
 */

#ifndef ANALYSEYIELD
#define ANALYSEYIELD

#include<map>
#include<string>
#include"BinVector.h"
#include"Analyse.h"

class AnalyseYield: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     */
    AnalyseYield(TreeWrapper *Tree);
    /**
     * Get the background subtracted yield
     */
    double GetBackgroundSubtractedYield(int Bin) const;
    /**
     * Run the analysis and save the double tag yields
     * @param Filename Filename where the yields are saved
     */
    void CalculateDoubleTagYields(const std::string &Filename);
    /**
     * Function that saves the results to a file
     * @param Filename Filename where the yields are saved
     */
    void SaveResults(const std::string &Filename) const;
  private:
    /**
     * Map with the yields in each bin and each region
     */
    std::map<char, BinVector<double>> m_Yields;
    /**
     * Events outside of the MBC region
     */
    int m_EventsOutsideMBCSpace;
    /**
     * Events outside of phase space
     */
    int m_EventsOutsidePhaseSpace;
};

#endif
