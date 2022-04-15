#ifndef CROPPED_VOR_DIAG_H
#define CROPPED_VOR_DIAG_H

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Identity_policy_2.h>
#include <CGAL/Polygon_2.h>

//#include "/Users/siviigni/Code/Balls/UnionOfBalls/src/ConvConvIntersection.h"

// put VoronoiDiagram and Polygon type as templates
// and check that the derived type Point_2 is the same
template <typename TKernel, typename TGraph>  
  class CroppedVoronoiDiagram
{
  typedef TKernel Kernel;
  typedef TGraph Graph;
  
  typedef CGAL::Delaunay_triangulation_2<Kernel> DT;
  typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
  typedef CGAL::Identity_policy_2<DT,AT> AP; // adaptation policy --> here keep everything
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP> VoronoiDiagram;

  typedef CGAL::Polygon_2<Kernel> Polygon_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef CGAL::Segment_2<Kernel> Segment_2;
  typedef CGAL::Ray_2<Kernel> Ray_2;
  typedef CGAL::Line_2<Kernel> Line_2;
  typedef CGAL::Iso_rectangle_2<Kernel> Iso_rectangle_2;
  
  /* // Graph */
  
  /* struct VertexProperty  */
  /* {  */
  /*   Point_2 id; */
  /* }; */
  
  
  /* typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexProperty> Graph; */
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex_t; 
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_t; 
  
  
  // ----------------------- Standard services ------------------------------
 public:
  
  /**
   * Destructor.
   */
  ~CroppedVoronoiDiagram();
  
  
  /**
   * Default Constructor.
   */
  CroppedVoronoiDiagram(const VoronoiDiagram & VD, const Polygon_2 &P);
  
  /**
   * Copy constructor.
   * @param other the object to clone.
   */
  CroppedVoronoiDiagram( const CroppedVoronoiDiagram & other );
  
  
    /**
   * Assignment.
   * @param other the object to copy.
   * @return a reference on 'this'.
   */
  CroppedVoronoiDiagram & operator=( const CroppedVoronoiDiagram & other );
  
  //void IntersectionSegmentPolygon(const Segment_2 & seg, const Polygon_2 & P, std::vector<Point_2> & inter);
  void intersectionSegmentPolygon(const Segment_2 & seg, const Polygon_2 & P, std::vector<Point_2> & inter);
  
  void cropSegment(const Segment_2 & seg, const Polygon_2 & P, std::vector<Segment_2> & cropped);
  bool computeSegFromHalfedge(const typename VoronoiDiagram::Halfedge & e, const Iso_rectangle_2 & R, Segment_2 & seg);
  
  void computeCroppedVoronoiDiagram();
  void computeGraph();

  void computeNeighbours();
  void initMap();
  
  std::vector<Segment_2> getCroppedVoronoiDiagram();
  Graph getVDGraph();

  
  void drawVDGraph(DGtal::Board2D & board);
  
  
 protected :
    VoronoiDiagram myVD;
    Polygon_2 myP;
    std::vector<std::pair<Segment_2,Point_2> > myCroppedVD2;
    
    std::vector<Segment_2> myCroppedVD;

    std::vector<std::vector<Segment_2>> croppedFaces;
    
    Graph myVDGraph; 
    
 private :
  std::map<Point_2,int> myMap; // map from Voronoi sites to inices
  
  std::vector<std::vector<int>> neighbours; // list of neighbour sites indices for each site
  
  

    std::map<Point_2,Vertex_t> mapGraphVertices; // map from Voronoi vertices coordinates to graph vertices
    
    
};



///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods.
#include "CroppedVoronoiDiagram.ih"


#endif //CROPPED_VOR_DIAG_H
