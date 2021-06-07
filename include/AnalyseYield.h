// Martin Duy Tat 19th May 2021
/**
 * AnalyseYield is a class for calculating the double tag yield from reconstructed events
 */

#ifndef ANALYSEYIELD
#define ANALYSEYIELD

#include<map>
#include<utility>
#include<vector>
#include<string>
#include"BinVector.h"
#include"Analyse.h"
#include"TMatrixT.h"

class AnalyseYield: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     * @param SubtractBackground Set to true to subtract off sideband and peaking backgrounds
     * @param Optional argument with file name to specify peaking backgrounds in each bin
     */
    AnalyseYield(TreeWrapper *Tree, bool SubtractBackground = true, const std::string &PeakingBackgroundFile = std::string(""));
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
    /**
     * Function that plots the Dalitz distribution of each bin and saves it as a .png file
     * @param Filename First part of the filename of the .png files, the bin number of file extension is appended
     */
    void SaveDalitzDistributions(const std::string &Filename) const;
  private:
    /**
     * Map with the yields in each bin and each region
     */
    std::map<char, BinVector<double>> m_Yields;
    /**
     * Vector of Dalitz positions of all events in each bin
     */
    BinVector<std::vector<std::pair<double, double>>> m_DalitzCoordinates;
    /**
     * If true, sideband subtraction is used to remove backgrounds
     */
    bool m_SubtractBackground;
    /**
     * List of peaking background in each bin and each region
     */
    std::map<char, BinVector<double>> m_PeakingBackground;
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
