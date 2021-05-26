// Martin Duy Tat 26th May 2021
/**
 * AnalysePeakingBackgrounds is a class for counting the binned yield of peaking backgrounds in inclusive MC
 */

#ifndef ANALYSEPEAKINGBACKGROUNDS
#define ANALYSEPEAKINGBACKGROUNDS

#include<map>
#include<string>
#include"BinVector.h"
#include"Analyse.h"
#include"TreeWrapper.h"

class AnalysePeakingBackgrounds: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     * @param Filename of text file with a list of iDcyTr numbers labelling the peaking backgrounds
     */
    AnalysePeakingBackgrounds(TreeWrapper *Tree, const std::string &Filename);
    /**
     * Function that count the number of peaking background events in each bin
     * @param Filename of file where peaking backgrounds are saved
     */
    void CalculatePeakingBackgrounds(const std::string &Filename);
    /**
     * Function that saves the list of peaking backgrounds to a file
     * @param Filename of file where peaking backgrounds are saved
     */
    void SavePeakingBackgrounds(const std::string &Filename) const;
  private:
    /**
     * Map containing all iDcyTr numbers and the binned yields
     */
    std::map<int, BinVector<int>> m_PeakingBackgrounds;
};

#endif
