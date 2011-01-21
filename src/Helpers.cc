#include "../include/Helpers.h"
#include <cstdlib>
#include <stdexcept>

namespace skylens {

  // general helper functions

  std::string getDatapath() {
    char* cp = getenv("SKYLENSDATAPATH");
    if (cp == NULL)
      throw std::runtime_error("SkyLens: environment variable SKYLENSDATAPATH is not set!");
    return std::string(cp);
  }

  void test_open(std::ifstream& ifs, std::string datapath, std::string& filename) {
    ifs.close();
    ifs.open((datapath+"/"+filename).c_str());
    if (ifs.fail()) {
      ifs.close();
      ifs.open(filename.c_str());
      if (ifs.fail()) 
	throw std::runtime_error("SkyLens: " + filename + " missing!");
    } else
      filename = datapath+"/"+filename;
  }

  std::vector<std::string> split(std::string s, char delimiter) {
    std::vector<std::string> chunks;
    std::string::size_type front = s.find(delimiter,0);
    if (front != std::string::npos) {
      std::string::size_type back = 0;
      while (front != std::string::npos) {
	chunks.push_back(s.substr(back,front-back));
	back = front+1;
	front = s.find(delimiter,back);
      }
      front = s.size();
      chunks.push_back(s.substr(back,front-back));
    } else 
      chunks.push_back(s);
    return chunks;
  }


} // end namespace
