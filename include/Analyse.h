// Martin Duy Tat 19th May 2021
/**
 * Analyse is a base class for running analyses
 */

#ifndef ANALYSE
#define ANALYSE

#include"TreeWrapper.h"
#include"BinningScheme.h"

class Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     */
    Analyse(TreeWrapper *Tree);
  protected:
    /**
     * Convert array index
     */
    int ArrayIndex(int i) const;
    /**
     * Determine MBC region for the current event
     */
    char DetermineMBCRegion() const;
    /**
     * Get reconstructed bin number for the current event
     */
    int DetermineReconstructedBinNumber() const;
    /**
     * Get generated bin number
     */
    int DetermineGeneratorBinNumber() const;
    /**
     * Tree with data
     */
    TreeWrapper *m_Tree;
    /**
     * Binning scheme
     */
    BinningScheme m_BinningScheme;
};

#endif
