// Martin Duy Tat 18th May 2021

#include<stdexcept>
#include<string>
#include<fstream>
#include<sstream>
#include"Settings.h"

KpiSettings& KpiSettings::Get() {
  static KpiSettings instance;
  return instance;
}

void KpiSettings::Initialize(const std::string &Filename) {
  std::ifstream SettingsFile(Filename);
  std::string line;
  while(std::getline(SettingsFile, line)) {
    if(line[0] == '#') {
      continue;
    }
    std::string Setting, Value;
    std::stringstream ss(line);
    ss >> Setting >> Value;
    m_Settings.insert({Setting, Value});
  }
}

double KpiSettings::GetDouble(const std::string &Setting) const {
  auto Search = m_Settings.find(Setting);
  if(Search == m_Settings.end()) {
    throw std::out_of_range(std::string("Could not find the setting ") + Setting);
  }
  return std::stod(m_Settings.at(Setting));
}

int KpiSettings::GetInt(const std::string &Setting) const {
  auto Search = m_Settings.find(Setting);
  if(Search == m_Settings.end()) {
    throw std::out_of_range(std::string("Could not find the setting ") + Setting);
  }
  return std::stoi(m_Settings.at(Setting));
}

std::string KpiSettings::GetString(const std::string &Setting) const {
  auto Search = m_Settings.find(Setting);
  if(Search == m_Settings.end()) {
    throw std::out_of_range(std::string("Could not find the setting ") + Setting);
  }
  return m_Settings.at(Setting);
}

bool KpiSettings::GetBool(const std::string &Setting) const {
  auto Search = m_Settings.find(Setting);
  if(Search == m_Settings.end()) {
    throw std::out_of_range(std::string("Could not find the setting ") + Setting);
  }
  return m_Settings.at(Setting) == "true";
}

KpiSettings::KpiSettings() {
}
