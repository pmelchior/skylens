#ifndef SKYLENS_HELPERS_H
#define SKYLENS_HELPERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace skylens {

#define TOSTRING(s) #s
#define STRINGIFY(s) TOSTRING(s)

  /// Get path to data files.
  /// This is taken from the environment variable \p SKYLENSDATAPATH.
  std::string getDatapath();
  /// Test and open a file.
  /// Tries to open the file \p filename either located in the 
  /// \p datapath directory or in "/". On success, \p ifs is a valid
  /// \p std::ifstream and \p filename  will be set to
  /// the name of the opened file. If the file cannot be found/opened, an
  /// exception of type \p std::runtime_error is thrown.
  void test_open(std::ifstream& ifs, std::string datapath, std::string& filename);
  /// Split string \p s into substrings, separated by \p delimiter.
  std::vector<std::string> split(std::string s, char delimiter);
}

#endif
