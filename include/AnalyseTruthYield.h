// Martin Duy Tat 19th May 2021
/**
 * AnalyseTruthYield goes through the MC truth tuple and count the number of events generated in each bin
 */

#ifndef ANALYSETRUTHYIELD
#define ANALYSETRUTHYIELD

#include<string>
#include"Analyse.h"
#include"BinVector.h"
#include"TreeWrapper.h"

class AnalyseTruthYield: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     */
    AnalyseTruthYield(TreeWrapper *Tree);
    /**
     * Run the analysis and save true double tag efficiencies
     * @param Filename Filename of text file to save double tag efficiencies
     */
    void CalculateTruthYield(const std::string &Filename);
  private:
    /**
     * Double tag yields from the generator
     */
    BinVector<double> m_GeneratorYields;
    /**
     * Number of events outside of phase space
     */
    int m_EventsOutsidePhaseSpace;
};

#endif
