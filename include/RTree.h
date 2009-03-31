#ifndef SKYLENS_RTREE_H
#define SKYLENS_RTREE_H

#include <spatialindex/SpatialIndex.h>
#include <frame/Catalog.h>
#include <frame/Point2D.h>
#include <vector>
#include <list>

namespace skylens {
  class RTree {
  public:
    RTree();
    RTree(const std::vector<SpatialIndex::Region>& vr);
    ~RTree();
    const std::list<unsigned int>& getMatches(const shapelens::Point2D<double>& P);

  private:
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
      const std::list<unsigned int>& getList() const  {
	return l;
      }
    private:
      std::list<unsigned int> l;
    };
    ListVisitor lvis;
    SpatialIndex::ISpatialIndex* tree;
    SpatialIndex::IStorageManager* mem;
  };
} // end namespace

#endif
