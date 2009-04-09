#include <RTree.h>

using namespace shapelens;
using namespace skylens;

RTree::RTree(const std::vector<shapelens::Rectangle<double> >& vr) {
  mem = SpatialIndex::StorageManager::createNewMemoryStorageManager();
  SpatialIndex::id_type indexIdentifier;
  tree = SpatialIndex::RTree::createNewRTree(*mem, 0.7, 20, 20, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
  // add nodes to tree
  for (unsigned long i=0; i < vr.size(); i++) {
    SpatialIndex::Region r(vr[i].ll.c_array(),vr[i].tr.c_array(), 2);
    tree->insertData(0, 0, r, i);
  }
}

RTree::~RTree() {
  delete tree;
  delete mem;
}

const std::list<unsigned long>& RTree::getMatches(const shapelens::Point2D<double>& P) {
  lvis.clear();
  SpatialIndex::Point test(P.c_array(),2);
  tree->pointLocationQuery(test,lvis);
  return lvis.getList();
}
