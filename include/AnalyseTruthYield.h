// Martin Duy Tat 19th May 2021
/**
 * AnalyseTruthYield goes through the MC truth tuple and count the number of events generated in each bin
 */

#ifndef ANALYSETRUTHYIELD
#define ANALYSETRUTHYIELD

#include<string>
#include"BinVector.h"
#include"TreeWrapper.h"

class AnalyseTruthYield {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     */
    AnalyseTruthYield(TreeWrapper *Tree);
    /**
     * Run the analysis and save true double tag efficiencies
     * @param SignalMCYield The yield from signal MC
     * @param Filename Filename of text file to save double tag efficiencies
     * @return Returns the bin efficiencies
     */
    BinVector<double> CalculateTruthYield(const BinVector<double> &SignalMCYield, const std::string &Filename);
  private:
    /**
     * Double tag yields from the generator
     */
    BinVector<double> m_GeneratorYields;
};

#endif
