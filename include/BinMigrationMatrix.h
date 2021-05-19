// Martin Duy Tat 19th May 2021
/**
 * BinMigrationMatrix is a class representing the information about the bin migration of double tags, and it can be used to correct raw yields
 */

#ifndef BINMIGRATIONMATRIX
#define BINMIGRATIONMATRIX

#include"TMatrix.h"

class BinMigrationMatrix {
  public:
    /**
     * Constructor that sets up the matrix
     * @param NBins Number of bins in binning scheme (not counting negative bins)
     */
    BinMigrationMatrix(int NBins);
    /**
     * Function that takes care of indexing
     * @param BinNumber The bin number
     * @return The array index (negative indices are places at the end)
     */
    int BinIndex(int BinNumber) const;
    /**
     * Function for adding an event based on its original and reconstructed bin number
     */
    void AddBinnedEvent(int ReconstructedBinNumber, int GeneratorBinNumber);
    /**
     * Normalize the rows of the bin migration matrix
     */
    TMatrixT<double> Normalize() const;
    /**
     * Correct yields for bin migration
     * @param OriginalYields BinVector with raw yields
     * @return Yields corrected for bin migration
     */
    BinVector<double> CorrectYields(const BinVector<double> &OriginalYields) const;
  private:
    /**
     * Bin migration matrix
     */
    TMatrixT<double> m_BinMigration;
};

#endif
