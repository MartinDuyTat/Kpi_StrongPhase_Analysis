// Martin Duy Tat 21st July 2021
/** 
 * PresentYields is a class that plots the Kpi yields and the prediction for comparison
 */

#ifndef PRESENTYIELDS
#define PRESENTYIELDS

#include<string>
#include"DoubleTagYields.h"
#include"HadronicParameters.h"

class PresentYields {
  public:
    /**
     * Constructor that stores all the relevant data
     * @param K0pipiKpiYieldFile File with K0pipi vs Kpi yields
     * @param K0KKKpiYieldFile File with K0KK vs Kpi yields
     * @param K0pipicisiHadronicFile File with K0pipi strong phase \f$c_i, s_i\f$ parameters
     * @param K0pipicisiHadronicFile File with K0pipi hadronic \f$K_i\f$ parameters
     * @param K0KKHadronicFile File with K0KK strong phase \f$c_i, s_i\f$ parameters
     * @param K0KKHadronicFile File with K0KK hadronic \f$K_i\f$ parameters
     * @param K0Mode "KS" or "KL"
     * @param rDcosDelta Assumed value of \f$r_D^{K\pi}\cos(\delta_D^{K\pi})\f$
     * @param rDsinDelta Assumed value of \f$r_D^{K\pi}\sin(\delta_D^{K\pi})\f$
     * @param DrawK0KK Set to true to also draw K0KK yields
     */
    PresentYields(const std::string &K0pipiKpiYieldFile, const std::string &K0KKKpiYieldFile, const std::string &K0pipicisiHadronicFile, const std::string &K0pipiKiHadronicFile, const std::string &K0KKcisiHadronicFile, const std::string &K0KKKiHadronicFile, const std::string &K0Mode, double rDcosDelta, double rDsinDelta, bool DrawK0KK = true);
    /**
     * Function that plots the Kpi yield and the prediction
     * @param Filename Filename to save plot to
     */
    void PlotYieldPresentation(const std::string &Filename) const;
  private:
    /**
     * K0pipi number of bins
     */
    int m_K0pipiBins;
    /**
     * K0KK number of bins
     */
    int m_K0KKBins;
    /**
     * The K0pipi vs Kpi double tag yield data
     */
    DoubleTagYields m_K0pipiDoubleTagYield;
    /**
     * The K0KK vs Kpi double tag yield data
     */
    DoubleTagYields m_K0KKDoubleTagYield;
    /**
     * The K0pipi hadronic parameters
     */
    HadronicParameters m_K0pipiHadronicParameters;
    /**
     * The K0KK hadronic parameters
     */
    HadronicParameters m_K0KKHadronicParameters;
    /**
     * "KS" or "KL"
     */
    std::string m_K0Mode;
    /**
     * \f$r_D^{K\pi}\cos(\delta_D^{K\pi})\f$
     */
    double m_rDcosDelta;
    /**
     * \f$r_D^{K\pi}\sin(\delta_D^{K\pi})\f$
     */
    double m_rDsinDelta;
};

#endif
