#ifndef SKYLENS_H
#define SKYLENS_H

#include "Helpers.h"
#include "Conversion.h"
#include "Layer.h"
#include "LensingInformation.h"
#include "Observation.h"
#include "PSF.h"
#include "RNG.h"
#include "RTree.h"
#include "SourceCatalog.h"
#include "Telescope.h"
#include "Filter.h"
#include "SED.h"
#include "Cosmology.h"

namespace skylens {
  /// Singleton'ed access the application DB.
  typedef Singleton<SQLiteDB> ApplicationDB;
}

#endif
