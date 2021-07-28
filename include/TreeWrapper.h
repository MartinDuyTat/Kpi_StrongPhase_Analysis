// Martin Duy Tat 19th May 2021
/**
 * TreeWrapper is a wrapper class containing the NTuple with double tags and its information is more easily accessed
 */

#ifndef TREEWRAPPER
#define TREEWRAPPER

#include<vector>
#include<string>
#include"TChain.h"
#include"TEntryList.h"
#include"TLorentzVector.h"

/**
 * ReconstructedKinematics is a simple struct containing all variables describing the reconstructed kinematics
 */
struct ReconstructedKinematics {
  /**
   * Flag that is 1 is the Kalman kinematic fit is a success
   */
  int KalmanFitSuccess;
  /**
   * Charge of tag kaon
   */
  int TagKCharge;
  /**
   * Beam constrained mass on the signal side
   */
  double SignalMBC;
  /**
   * Missing mass squared on the signal side
   */
  double SignalMMiss2;
  /**
   * Beam constrained mass on the signal side
   */
  double TagMBC;
  /**
   * Missing energy minus missing momentum on the tag side
   */
  double TagUMiss;
  /**
   * \f$K^0\f$ momentum
   */
  TLorentzVector K0_P;
  /**
   * \f$K^0\f$ momentum after Kalman kinematic fit
   */
  TLorentzVector K0Kalman_P;
  /**
   * \f$K^+\f$ momentum
   */
  TLorentzVector KPlus_P;
  /**
   * \f$K^+\f$ momentum after Kalman kinematic fit
   */
  TLorentzVector KPlusKalman_P;
  /**
   * \f$K^-\f$ momentum
   */
  TLorentzVector KMinus_P;
  /**
   * \f$K^-\f$ momentum after Kalman kinematic fit
   */
  TLorentzVector KMinusKalman_P;
};

/**
 * GeneratorKinematics is a simple struct containing all variables describing the generator kinematics
 */
struct GeneratorKinematics {
  /**
   * Number of particles generated in event
   */
  int NumberParticles;
  /**
   * List of particle IDs
   */
  std::vector<int> ParticleIDs = std::vector<int>(100);
  /**
   * List of mother indices
   */
  std::vector<int> MotherIndex = std::vector<int>(100);
  /**
   * List of generated momenta in the x-direction
   */
  std::vector<double> TruePx = std::vector<double>(100);
  /**
   * List of generated momenta in the y-direction
   */
  std::vector<double> TruePy = std::vector<double>(100);
  /**
   * List of generated momenta in the z-direction
   */
  std::vector<double> TruePz = std::vector<double>(100);
  /**
   * List of generated energies in the x-direction
   */
  std::vector<double> TrueEnergy = std::vector<double>(100);
  /**
   * Decay topology number from TopoAna
   */
  int iDcyTr = -1;
};

class TreeWrapper {
  public:
    /**
     * Constructor that sets up the NTuple
     * @param Filename Name of text file containing a list of ROOT file(s) (can be multiple files using wild card)
     * @param TreeName Name of ROOT object
     * @param CutFile File containing selection cuts
     * @param DataType "Data", "SignalMC" or "TruthTuple"
     * @param MomentumSmearing Smear the K0KK momenta by this amount
     */
    TreeWrapper(const std::string &Filename, const std::string &TreeName, const std::string &CutFile, const std::string &DataType, double MomentumSmearing = 0.0);
    /**
     * Function that sets all the branch addresses
     * @param DataType "Data", "SignalMC" or "TruthTuple"
     */
    void SetBranchAddresses(const std::string &DataType);
    /**
     * Get reconstructed kinematics
     */
    const ReconstructedKinematics& GetReconstructedKinematics() const;
    /**
     * Get generator kinematics
     */
    const GeneratorKinematics& GetGeneratorKinematics() const;
    /**
     * Get number of events after cuts
     */
    int GetEntries() const;
    /**
     * Get entry
     * @param i Entry number
     */
    void GetEntry(int i);
    /**
     * Get signal mode
     */
    std::string GetSignalMode() const;
    /**
     * Get tag mode
     */
    std::string GetTagMode() const;
    /**
     * Smear the K0KK momenta
     */
    void SmearMomenta();
  private:
    /**
     * TChain object containing all data
     */
    TChain m_Chain;
    /**
     * TEntryList object that only selects events passing the selection cuts
     */
    TEntryList *m_elist;
    /**
     * Signal mode, "KSKK" or "KLKK"
     */
    std::string m_SignalMode;
    /**
     * Tag mode, "Kpi", "Kpipi0", "Kpipipi" or "KeNu"
     */
    std::string m_TagMode;
    /**
     * Struct with reconstructed kinematics
     */
    ReconstructedKinematics m_RecKinematics;
    /**
     * Struct with generator kinematics
     */
    GeneratorKinematics m_GenKinematics;
    /**
     * Smear the K0KK momenta by a random Gaussian with this standard deviation, set to zero for no smearing
     */
    double m_MomentumSmearing;
};

#endif
