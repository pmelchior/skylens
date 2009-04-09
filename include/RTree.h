#ifndef SKYLENS_RTREE_H
#define SKYLENS_RTREE_H

#include <spatialindex/SpatialIndex.h>
#include <frame/Catalog.h>
#include <frame/Point2D.h>
#include <modelfit/SourceModel.h>
#include <vector>
#include <list>

namespace skylens {
  /// Class for efficient 2D lookups.
  /// Uses a R* tree index from http://trac.gispython.org/spatialindex/wiki/Releases.
  class RTree {
  public:
    /// Default constructor.
    RTree();
    /// Constructor from a set of Rectangle entities.
    RTree(const std::vector<shapelens::Rectangle<double> >& vr);
    /// Destructor.
    ~RTree();
    /// Find object whose support Rectangle overlaps with \p P.
    /// The list contains the vector indices of those objects whose rectangles
    /// were given at construction time.
    const std::list<unsigned long>& getMatches(const shapelens::Point2D<double>& P);

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
      const std::list<unsigned long>& getList() const  {
	return l;
      }
    private:
      std::list<unsigned long> l;
    };
    ListVisitor lvis;
    SpatialIndex::ISpatialIndex* tree;
    SpatialIndex::IStorageManager* mem;
  };
} // end namespace

#endif
