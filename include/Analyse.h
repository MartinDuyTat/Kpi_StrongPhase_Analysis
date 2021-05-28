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
     * Determine energy region for the current event
     * For fully reconstructed double tags, this is the 2D MBC plane
     * For partially reconstructed tags, the reconstructed tag must have MBC between 1.86 and 1.87 GeV, and the missing mass/energy for the partially reconstructed tag
     */
    char DetermineEnergyRegion() const;
    /**
     * Get reconstructed bin number for the current event
     */
    int DetermineReconstructedBinNumber() const;
    /**
     * Determine closest bin number, assuming the event is outside the Dalitz phase space
     */
    int DetermineMappedReconstructedBinNumber() const;
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
