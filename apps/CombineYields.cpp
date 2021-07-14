// Martin Duy Tat 15th July 2021
/**
 * CombineYields is an application that combines the yields and their statistical and systematic errors
 * @param 1 Filename of combined yields output
 * @param 2, 3, ... Filenames of text files with yields to be combined
 */

#include<iostream>
#include<vector>
#include<string>
#include"YieldCombiner.h"

int main(int argc, char *argv[]) {
  if(argc < 3) {
    std::cout << "Need more than 2 input arguments\n";
    return 0;
  }
  std::cout << "Combining yields from " << argc - 2 << " dataset(s)\n";
  std::vector<std::string> Filenames;
  for(int i = 2; i < argc; i++) {
    std::string Filename(argv[i]);
    std::cout << "Adding yields from " << Filename << "\n";
    Filenames.push_back(Filename);
  }
  std::cout << "Combining yields...\n";
  YieldCombiner Combiner(2, Filenames);
  Combiner.SaveYields(std::string(argv[1]));
  std::cout << "Final yields ready!\n";
  return 0;
}
