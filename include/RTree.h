#ifndef SKYLENS_RTREE_H
#define SKYLENS_RTREE_H

#include <spatialindex/SpatialIndex.h>
#include <shapelens/frame/Catalog.h>
#include <shapelens/frame/Point.h>
#include <shapelens/modelfit/SourceModel.h>
#include <vector>
#include <list>

namespace skylens {
  /// Class for efficient 2D lookups.
  /// Uses a R* tree index from http://trac.gispython.org/spatialindex/wiki/Releases.
  class RTree {
  public:
    /// Default constructor.
    RTree();
    /// Insert nodes.
    /// Nodes are Rectangle entities for which the RTree stores their index.
    void insertNodes(const std::vector<shapelens::Rectangle<double> >& vr);
    /// Destructor.
    ~RTree();
    /// Find object whose support Rectangle overlaps with \p P.
    /// The list contains the vector indices of those objects whose rectangles
    /// were given at construction time.
    std::list<unsigned long> getMatches(const shapelens::Point<double>& P) const;

  private:
    SpatialIndex::ISpatialIndex* tree;
    SpatialIndex::IStorageManager* mem;
  };
} // end namespace

#endif
