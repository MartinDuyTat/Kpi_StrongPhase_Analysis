// Martin Duy Tat 13th May 2021

#ifndef DOUBLETAGYIELDS
#define DOUBLETAGYIELDS

#include<string>
#include"BinVector.h"

/**
 * DoubleTagYields is a class that stores the double tag yields of each bin and its error
 */
class DoubleTagYields {
  public:
    /**
     * Constructor that reads the yields from a file
     * @param Filename Filename of text files with double tag yields
     * @param NBins Number of bins
     */
    DoubleTagYields(const std::string &Filename, int NBins);
    /**
     * Get yield
     * @param Bin Bin number
     */
    double GetYield(int Bin) const;
    /**
     * Get yield error
     * @param Bin Bin number
     */
    double GetYieldError(int Bin) const;
    /**
     * Function that prints the yields
     */
    void PrintYields() const;
  private:
    /**
     * Vector with yields in each bin
     */
    BinVector<double> m_Yield;
    /**
     * Vector with yield errors in each bin
     */
    BinVector<double> m_YieldError;
};

#endif
