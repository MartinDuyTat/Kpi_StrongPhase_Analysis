// Martin Duy Tat 18th May 2021

#ifndef KPISETTINGS
#define KPISETTINGS

#include<map>
#include<string>

class KpiSettings {
  public:
    /**
     * Get settings object
     */
    static KpiSettings& Get();
    /**
     * Initialize settings with options from text file
     * @param Filename Filename of text file with settings
     */
    void Initialize(const std::string &Filename);
    /**
     * Get options in double format
     * @param Setting Name of the setting
     */
    double GetDouble(const std::string &Setting);
    /**
     * Get options in int format
     * @param Setting Name of the setting
     */
    int GetInt(const std::string &Setting);
    /**
     * Get options in string format
     * @param Setting Name of the setting
     */
    std::string GetString(const std::string &Setting);
  private:
    /**
     * Private constructor because I only want one copy of this object around
     */
    KpiSettings();
    /**
     * Delete copy constructor
     */
    KpiSettings(const KpiSettings &Settings) = delete;
    /**
     * Delete assignment operator
     */
    KpiSettings& operator=(const KpiSettings &Settings) = delete;
    /**
     * Map that contains the settings name and value, both in string format
     */
    std::map<std::string, std::string> m_Settings;
};

#endif
