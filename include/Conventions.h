#ifndef SKYLENS_CONVENTIONS_H
#define SKYLENS_CONVENTIONS_H

#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)
#include <string>
namespace skylens {
  const std::string datapath = TOSTRING(DATAPATH);
}

#endif
