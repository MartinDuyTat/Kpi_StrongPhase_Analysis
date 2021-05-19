// Martin Duy Tat 19th May 2021
/**
 * AnalyseBinMigration is a class for calculating the bin migration matrix using signal MC
 */

#ifndef ANALYSEBINMIGRATION
#define ANALYSEBINMIGRATION

#include<string>
#include"Analyse.h"
#include"TreeWrapper.h"
#include"TMatrixT.h"

class AnalyseBinMigration: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     */
    AnalyseBinMigration(TreeWrapper *Tree);
    /**
     * Run the analysis and get the double tag yields in each reconstructed and generated bin
     * @param Filename Filename where bin yields are saved in matrix form
     */
    void CalculateBinMigrationYields(const std::string &Filename);
    /**
     * Function that saves the results to a file
     * @param Filename Filename where bin yields are saved in matrix form
     */
    void SaveResults(const std::string &Filename) const;
    /**
     * Get the bin migration matrix (normalize the yields and invert the matrix)
     */
    TMatrixT<double> GetBinMigrationMatrix() const;
  private:
    /**
     * Matrix of generated bin number (rows) and reconstructed bin number (column) with positive bins first then negative bins (+1 +2 -1 -2)
     */
    TMatrixT<double> m_BinYields;
};

#endif
