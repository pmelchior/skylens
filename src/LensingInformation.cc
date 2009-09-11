#include "../include/LensingInformation.h"
#include <limits>

namespace skylens {
  LensingInformation::LensingInformation() :
    z_first_lens (std::numeric_limits<double>::infinity()) {}
}
