// Martin Duy Tat 26th May 2021
/**
 * AnalysePeakingBackgrounds is a class for counting the binned yield of peaking backgrounds in inclusive MC
 * The iDcyTr numbers from TopoAna are read from a file, the first line contains the signal components, which are ignored, the second line has all the peaking backgrounds and any numbers not listed are put in some "Other" category
 * The third line contains the scaling factors, other than the standard luminosity scale factor, which is read from the settings file, and the last entry is for the "Other" backgrounds category
 */

#ifndef ANALYSEPEAKINGBACKGROUNDS
#define ANALYSEPEAKINGBACKGROUNDS

#include<map>
#include<utility>
#include<string>
#include<vector>
#include"BinVector.h"
#include"Analyse.h"
#include"TreeWrapper.h"

class AnalysePeakingBackgrounds: public Analyse {
  public:
    /**
     * Constructor that saves a pointer to the NTuple containing the data and also sets up binning scheme
     * @param Tree Pointer to data
     * @param Filename of text file with a list of iDcyTr numbers labelling the signal and peaking backgrounds
     */
    AnalysePeakingBackgrounds(TreeWrapper *Tree, const std::string &Filename);
    /**
     * Function that count the number of peaking background events in each bin
     * @param Filename of file where peaking backgrounds are saved
     * @param MCScale Luminosity scale of MC sample
     */
    void CalculatePeakingBackgrounds(const std::string &Filename, double MCScale = 1.0);
    /**
     * Function that saves the list of peaking backgrounds to a file
     * @param Filename of file where peaking backgrounds are saved
     */
    void SavePeakingBackgrounds(const std::string &Filename) const;
  private:
    /**
     * Map containing all iDcyTr numbers and the binned yields for each region
     */
    std::map<std::pair<int, char>, BinVector<double>> m_PeakingBackgrounds;
    /**
     * Scale factors to the yields of different modes to correct for different branching ratios in MC
     */
    std::map<int, double> m_ScaleFactors;
    /**
     * Binned yields of other backgrounds not listed
     */
    std::map<char, BinVector<double>> m_OtherBackgrounds;
    /**
     * List of signal iDcyTr numbers, these are ignored!
     */
    std::vector<int> m_SignalComponents;
};

#endif
