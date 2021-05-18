// Martin Duy Tat 13th May 2021

#ifndef BINVECTOR
#define BINVECTOR

#include<vector>
#include<stdexcept>
#include<fstream>
#include<sstream>

/**
 * BinVector is a wrapper class for std::vector but with negative bin indices and excluding zero
 */
template<class T>
class BinVector {
  public:
    /**
     * Default constructor
     */
    BinVector(): m_NegativeBins(0) {
  }
    /**
     * Constructor that allocates the correct number of vector elements
     * @param NegativeBins Set to true to include negative bin indices
     * @param Size
     */
    BinVector(bool NegativeBins, int Size): m_NegativeBins(NegativeBins), m_Size(Size), m_BinVector(std::vector<T>(NegativeBins ? 2*m_Size : m_Size)) {
    }
    /**
     * Operator overload of [] to account for negative bin indices
     * @param i Bin number
     */
    T& operator [](int i) {
      if(i == 0) {
	throw std::out_of_range("Bin number cannot be zero");
      }
      if(!m_NegativeBins && i < 0) {
	throw std::out_of_range("Bin number must be positive");
      }
      return m_BinVector[i > 0 ? i - 1 : -i - 1 + m_Size];
    }
    const T operator [](int i) const {
      if(i == 0) {
	throw std::out_of_range("Bin number cannot be zero");
      }
      if(!m_NegativeBins && i < 0) {
	throw std::out_of_range("Bin number must be positive");
      }
      return m_BinVector[i > 0 ? i - 1 : -i - 1 + m_Size];
    }
    /**
     * Function for getting size of binning
     */
    int Size() const {
      return m_Size;
    }
    /**
     * Iterator to the first element
     */
    typename std::vector<T>::iterator begin() {
      return m_BinVector.begin();
    }
    /**
     * Iterator to the last element
     */
    typename std::vector<T>::iterator end() {
      return m_BinVector.end();
    }
  private:
    /**
     * Flag that is true if negative bins are enabled
     */
    bool m_NegativeBins;
    /**
     * Number of bins
     */
    int m_Size;
    /**
     * Vector with data
     */
    std::vector<T> m_BinVector;
};

#endif
