#include "../include/RTree.h"

using namespace shapelens;
using namespace skylens;

// helper class
class ListVisitor : public SpatialIndex::IVisitor {
public:
  void visitNode(const SpatialIndex::INode& n) {}

  void visitData(const SpatialIndex::IData& d) {
    l.push_back(d.getIdentifier());
  }
  void visitData(std::vector<const SpatialIndex::IData*>& v) {
    for(std::vector<const SpatialIndex::IData*>::iterator iter = v.begin(); iter!= v.end(); iter++)
      l.push_back((*iter)->getIdentifier());
  }
  void clear() {
    l.clear();
  }
  const std::list<unsigned long>& getList() const  {
    return l;
  }
private:
  std::list<unsigned long> l;
};

RTree::RTree() {
  mem = SpatialIndex::StorageManager::createNewMemoryStorageManager();
  SpatialIndex::id_type indexIdentifier;
  tree = SpatialIndex::RTree::createNewRTree(*mem, 0.7, 20, 20, 2, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
}

void RTree::insertNodes(const std::vector<shapelens::Rectangle<double> >& vr) {
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

std::list<unsigned long> RTree::getMatches(const shapelens::Point2D<double>& P) const {
  ListVisitor lvis;
  SpatialIndex::Point test(P.c_array(),2);
  tree->pointLocationQuery(test,lvis);
  return lvis.getList();
}
