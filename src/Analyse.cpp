// Martin Duy Tat 19th May 2021

#include<stdexcept>
#include<vector>
#include<algorithm>
#include"Analyse.h"
#include"TreeWrapper.h"

Analyse::Analyse(TreeWrapper *Tree): m_Tree(Tree) {
}

int Analyse::ArrayIndex(int i) const {
  if(i > 0) {
    return i - 1;
  } else if (i < 0) {
    return -i + m_BinningScheme.GetNumberBins() - 1;
  } else {
    throw std::out_of_range("Cannot have bin number 0");
  }
}

char Analyse::DetermineMBCRegion() const {
  double SignalMBC = m_Tree->GetReconstructedKinematics().SignalMBC;
  double TagMBC = m_Tree->GetReconstructedKinematics().TagMBC;
  if(SignalMBC > 1.86 && SignalMBC < 1.87 && TagMBC > 1.86 && TagMBC < 1.87) {
    return 'S';
  } else if(SignalMBC > 1.86 && SignalMBC < 1.87 && TagMBC > 1.83 && TagMBC < 1.855) {
    return 'A';
  } else if(SignalMBC > 1.83 && SignalMBC < 1.855 && TagMBC > 1.86 && TagMBC < 1.87) {
    return 'B';
  } else if(SignalMBC > 1.83 && SignalMBC < 1.855 && TagMBC > 1.83 && TagMBC < 1.855) {
    if(TMath::Abs(SignalMBC - TagMBC) < 0.0035) {
      return 'C';
    } else if(TMath::Abs(SignalMBC - TagMBC) > 0.0055) {
      return 'D';
    } else {
      return 'F';
    }
  } else {
    return 'F';
  }
}

int Analyse::DetermineReconstructedBinNumber() const {
  ReconstructedKinematics RecKinematics = m_Tree->GetReconstructedKinematics();
  double M2Plus = (RecKinematics.KSKalman_P + RecKinematics.KPlusKalman_P).M2();
  double M2Minus = (RecKinematics.KSKalman_P + RecKinematics.KMinusKalman_P).M2();
  int KCharge = RecKinematics.TagKCharge;
  int BinNumber = m_BinningScheme.GetBinNumber(M2Plus, M2Minus, KCharge);
  if(BinNumber == 0) {
    M2Plus = (RecKinematics.KS_P + RecKinematics.KPlus_P).M2();
    M2Minus = (RecKinematics.KS_P + RecKinematics.KMinus_P).M2();
    return m_BinningScheme.GetBinNumber(M2Plus, M2Minus, KCharge);
  } else {
    return BinNumber;
  }
}

int Analyse::DetermineGeneratorBinNumber() const {
  std::vector<int> ParticleIDs = m_Tree->GetGeneratorKinematics().ParticleIDs;
  auto D0Position = std::find(ParticleIDs.begin(), ParticleIDs.end(), 421);
  auto D0barPosition = std::find(ParticleIDs.begin(), ParticleIDs.end(), -421);
  std::vector<int> D0Daughters = std::vector<int>(D0Position + 1, D0barPosition);
  std::vector<int> D0barDaughters = std::vector<int>(D0barPosition + 1, ParticleIDs.end());
  int KSKKMother;
  if(std::find(D0Daughters.begin(), D0Daughters.end(), 310) != D0Daughters.end()) {
    KSKKMother = 421;
  } else   if(std::find(D0barDaughters.begin(), D0barDaughters.end(), 310) != D0barDaughters.end()) {
    KSKKMother = -421;
  } else {
    throw std::logic_error("Cannot find KS daughter in ParticlesIDs");
  }
  int KSKKIndex = std::find(ParticleIDs.begin(), ParticleIDs.end(), KSKKMother) - ParticleIDs.begin();
  std::vector<double> TruePx = m_Tree->GetGeneratorKinematics().TruePx;
  std::vector<double> TruePy = m_Tree->GetGeneratorKinematics().TruePy;
  std::vector<double> TruePz = m_Tree->GetGeneratorKinematics().TruePz;
  std::vector<double> TrueEnergy = m_Tree->GetGeneratorKinematics().TrueEnergy;
  TLorentzVector KPlusP, KMinusP, KSP;
  int ParticlesFound = 0;
  for(int j = KSKKIndex + 1; j < m_Tree->GetGeneratorKinematics().NumberParticles; j++) {
    if(ParticleIDs[j] == -421) {
      break;
    }
    if(ParticleIDs[j] == 321) {
      KPlusP[0] = TruePx[j];
      KPlusP[1] = TruePy[j];
      KPlusP[2] = TruePz[j];
      KPlusP[3] = TrueEnergy[j];
      ParticlesFound++;
    } else if(ParticleIDs[j] == -321) {
      KMinusP[0] = TruePx[j];
      KMinusP[1] = TruePy[j];
      KMinusP[2] = TruePz[j];
      KMinusP[3] = TrueEnergy[j];
      ParticlesFound++;
    } else if(ParticleIDs[j] == 310) {
      KSP[0] = TruePx[j];
      KSP[1] = TruePy[j];
      KSP[2] = TruePz[j];
      KSP[3] = TrueEnergy[j];
      ParticlesFound++;
    }
  }
  if(ParticlesFound != 3) {
    throw std::logic_error("Could not find KSKK daughters");
  }
  double M2Plus = (KPlusP + KSP).M2();
  double M2Minus = (KMinusP + KSP).M2();
  int KCharge;
  if(KSKKMother == 421) {
    KCharge = std::find(D0barDaughters.begin(), D0barDaughters.end(), 321) != D0barDaughters.end() ? +1 : -1;
  } else {
    KCharge = std::find(D0Daughters.begin(), D0Daughters.end(), -321) != D0barDaughters.end() ? -1 : +1;
  }
  return m_BinningScheme.GetBinNumber(M2Plus, M2Minus, KCharge);
}
