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
#include"TMatrixT.h"

class AnalyseYield: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     * @param SubtractBackground Set to true to use a sideband subtraction
     */
    AnalyseYield(TreeWrapper *Tree, bool SubtractBackground = true);
    /**
     * Get the background subtracted yield
     */
    double GetBackgroundSubtractedYield(int Bin) const;
    /**
     * Run the analysis and save the double tag yields
     * @param BinMigrationMatrix Matrix that accounts for bin migration
     * @param Filename Filename where the yields are saved
     */
    void CalculateDoubleTagYields(const TMatrixT<double> &BinMigrationMatarix, const std::string &Filename);
    /**
     * Function that saves the results to a file
     * @param Filename Filename where the yields are saved
     */
    void SaveFinalYields(const TMatrixT<double> &BinMigrationMatrix, const std::string &Filename) const;
  private:
    /**
     * Map with the yields in each bin and each region
     */
    std::map<char, BinVector<double>> m_Yields;
    /**
     * If true, sideband subtraction is used to remove backgrounds
     */
    bool m_SubtractBackground;
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
