
#include <iostream>
#include <math.h>
#include <string>

// DGtal

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/DigitalSetBoundary.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"

// Boost Graph

#include <boost/graph/adjacency_list.hpp>

// CGAL

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>


#include "CroppedVoronoiDiagram.h"
#include "PowerMapZero.h"


using namespace DGtal;
using namespace std;


template <typename TKernel, typename DigitalSet>
DigitalSet digitize(const CGAL::Circle_2<TKernel> & b)
{
  typename TKernel::Point_2 c = b.center();
  typename TKernel::FT rr = b.squared_radius();
  
  double rr1 = rr.floatValue();
  Z2i::Point cc(c.x().floatValue(),c.y().floatValue());

  CGAL::Bbox_2 boudingbox = b.bbox(); 
  
  typename DigitalSet::Domain domain(Z2i::Point(floor(cc[0]-sqrt(rr1)),floor(cc[1]-sqrt(rr1))), Z2i::Point(ceil(cc[0]+sqrt(rr1)),ceil(cc[1]+sqrt(rr1))));
  
  DigitalSet digit(domain);
  
  for(auto const &p : domain)
    {
      typename TKernel::Point_2 pp(p[0],p[1]);
      if(b.has_on_positive_side(pp)) // open ball
	digit.insert(p);
    }
  return digit;
}


// Computes the circle corresponding to a vertex of the dag. 
template <typename  TKernel, typename DirectedGraph>
typename CGAL::Circle_2<TKernel> circleFromVertex(const typename DirectedGraph::vertex_descriptor & el, const DirectedGraph &dag)
{
  typename TKernel::Point_2 c(dag[el].id); // retrieve the center coordinates
  typename TKernel::Point_2 site(dag[el].site); // retrieve a closest site coordinates

  // compute the radius
  typename TKernel::FT squaredR = CGAL::squared_distance(c,site); 
   
  return CGAL::Circle_2<TKernel>(c,squaredR);
}
  


template <typename  TKernel, typename DirectedGraph>
void pencilFromVertex(const typename DirectedGraph::vertex_descriptor & el, const DirectedGraph &dag, typename CGAL::Circle_2<TKernel> & b1, typename CGAL::Circle_2<TKernel> &b2)
{
  b1 = circleFromVertex<TKernel,DirectedGraph>(el,dag);
  typename DirectedGraph::adjacency_iterator ibegin, iend;
  typename DirectedGraph::vertex_descriptor desc;
  for (boost::tie(ibegin, iend) = boost::adjacent_vertices(el, dag); ibegin != iend; ++ibegin)
    {
      desc = *ibegin;
    }
  if(desc == el)
    b2 = b1;
  else
    {
      b2 = circleFromVertex<TKernel,DirectedGraph>(desc,dag);
    }
}



template <typename  TKernel, typename DirectedGraph>
void pencilFromEdge(const typename DirectedGraph::edge_descriptor & el, const DirectedGraph &dag, typename CGAL::Circle_2<TKernel> & b1, typename CGAL::Circle_2<TKernel> &b2)
{
  typename DirectedGraph::vertex_descriptor vb1 = boost::source(el,dag);
  typename DirectedGraph::vertex_descriptor vb2 = boost::target(el,dag);

  b1 = circleFromVertex<TKernel,DirectedGraph>(vb1,dag);
  b2 = circleFromVertex<TKernel,DirectedGraph>(vb2,dag);

}

  

template <typename TKernel>
typename TKernel::FT power(const typename CGAL::Circle_2<TKernel> & b, const typename TKernel::Point_2 & x)
{
  typename TKernel::Point_2 c=b.center();
  typename TKernel::FT r= b.squared_radius();
  return (CGAL::squared_distance(c,x)-r);
}


// Given a point x and a pencil, compute the maximal ball of the pencil (lambda
// value) that contains x on its boundary. Assumes that x belongs to b1, not to b2
template <typename TKernel>
typename TKernel::FT computeLambdaOfRepBall(const typename TKernel::Point_2 &x, const typename CGAL::Circle_2<TKernel> & b1, const typename CGAL::Circle_2<TKernel> & b2)
{
  typename TKernel::FT lambda = power(b1,x)/(power(b1,x)-power(b2,x));
  return lambda;
}




// Given a lambda value and a pencil, compute the corresponding ball
// (center + radius)
template <typename TKernel>
typename CGAL::Circle_2<TKernel> ballOnPencil(const typename TKernel::FT & lambda, const typename CGAL::Circle_2<TKernel> & b1, const typename CGAL::Circle_2<TKernel> & b2)
{
  typename TKernel::Point_2 c1 = b1.center(), c2 = b2.center();
  typename TKernel::Point_2 c = CGAL::barycenter(c1,1-lambda,c2,lambda);
  
  typename TKernel::FT pow1 = power(b1, typename TKernel::Point_2(0,0));
  typename TKernel::FT pow2 = power(b2, typename TKernel::Point_2(0,0));
  
  
  typename TKernel::FT squaredR = CGAL::squared_distance(c,typename TKernel::Point_2(0,0)) - (1-lambda)*pow1 - lambda*pow2;

  return typename CGAL::Circle_2<TKernel>(c,squaredR); 
}



// Given a point x and a pencil, compute an antecedent of the maximal
// digital ball that contains x in its interior of the pencil. Assumes
// that x belongs to b1, not to b2  
template <typename TKernel, typename DigitalSet>
typename TKernel::FT computeLambdaOfCriticalBall(const typename TKernel::Point_2 &x, const typename CGAL::Circle_2<TKernel> & b1, const typename CGAL::Circle_2<TKernel> & b2)
{
  // compute Rep ball of x (supremum of all the balls that contain x)
  typename TKernel::FT lambda = computeLambdaOfRepBall<TKernel>(x,b1,b2);
  typename CGAL::Circle_2<TKernel> blambda = ballOnPencil<TKernel>(lambda,b1,b2);

   
  // compute the set of digital points in the ball
  DigitalSet digitBall = digitize<TKernel,DigitalSet>(blambda);
  DigitalSet digitBallb1 = digitize<TKernel,DigitalSet>(b1);
  
  typename TKernel::FT lambdaMax= 0; 

  // for all points in blambda but not in b1

  for(auto const &p : digitBall)
    {
      if(!digitBallb1(p))
	{
	  typename TKernel::Point_2 pp(p[0],p[1]);
	  typename TKernel::FT lambda1 = computeLambdaOfRepBall<TKernel>(pp,b1,b2);
	  if(lambda1 < 1 && lambda1 > lambdaMax)
	    lambdaMax = lambda1;
	}
    }

  if(lambdaMax != lambda)
    return (lambdaMax+lambda)/2;
  else
    return lambda/2;
}


template <typename Graph, typename DirectedGraph>
DirectedGraph directedGraphFromTree(const Graph &g, std::vector<typename DirectedGraph::vertex_descriptor> &topologicalOrdering, std::vector<typename DirectedGraph::edge_descriptor> &topologicalOrderingEdges)
{
  trace.beginBlock("Compute dag from tree");
  
  DirectedGraph dag;
  
  typedef typename Graph::vertex_descriptor Vertex_t_g;
  
  typedef typename DirectedGraph::vertex_descriptor Vertex_t_dag;
  typedef typename DirectedGraph::edge_descriptor Edge_t_dag;

  int n = num_vertices(g);
  
  std::queue<std::pair<Vertex_t_g,Vertex_t_dag> > queue;

  
  int r = rand() % n;
  typename Graph::vertex_iterator rootIt = vertices(g).first;
  for(int i=0;i<r;i++)
    {
      rootIt++;
    }
  
  Vertex_t_g root = *rootIt;
  
  Vertex_t_dag u_dag = add_vertex(typename DirectedGraph::vertex_property_type{g[root].id,g[root].site},dag);
  
  queue.push(std::make_pair(root,u_dag));

  
  topologicalOrdering.push_back(u_dag);
  
  std::cout << "root = " << g[root].id.x().floatValue() <<  " " << g[root].id.y().floatValue() << std::endl; 
  
  // Not efficient to find a point -> TODO : change for another data
  // structure with O(1) find cost --> sorted_vector with
  // std::lower_bound to find the element ??? 
  std::vector<Vertex_t_g> tagged;
  tagged.push_back(root);
  
  while(!queue.empty())
    {
      Vertex_t_g u_g = queue.front().first;
      Vertex_t_dag u_dag = queue.front().second;
      queue.pop();

      typename boost::graph_traits<Graph>::adjacency_iterator vi, vend;
      for(boost::tie(vi, vend) = adjacent_vertices(u_g,g); vi != vend; ++vi)
	{	  
	  if(std::find(tagged.begin(), tagged.end(), *vi) == tagged.end())
	    {
	      // add a vertex in the dag
	      Vertex_t_dag vi_dag = add_vertex(typename DirectedGraph::vertex_property_type{g[*vi].id,g[*vi].site},dag);
	      
	      // add an oriented edge in the dag between vi_dag and u_dag
	      	      
	      std::pair<Edge_t_dag, bool> e = add_edge(vi_dag,u_dag,dag);
	      	      
	      // add the pair (vi,vi_dag) to the queue
	      queue.push(std::make_pair(*vi,vi_dag));
	      tagged.push_back(*vi);

	      //set rank of vi_dag to k and decrement k
	      topologicalOrdering.push_back(vi_dag);
	      topologicalOrderingEdges.push_back(e.first);
	      
	    }
	}
      
    }
   
  trace.endBlock();
  return dag;
}

  

int main( int argc, char** argv )
{
  
  srand(time(0));

  int subsamp = 1;

  if(argc !=2 && argc != 4 )
    {
      std::cout << "Usage: " << argv[0] << " [input filename] (-sub d)\n" ;
      return 0;
    }
  else
    if(argc == 4)
      {
	if(strcmp(argv[2],"-sub")==0)
	  subsamp = atoi(argv[3]);
	else
	  {
	    std::cout << "Usage: " << argv[0] << " [input filename] (-sub d)\n" ;
	    return 0;
	  }
      }
  
  std::string imageFilename(argv[1]);

  
  //******** Part I: image -> Digital Set ********//

  // 1. Define Image type and read image file
  typedef DGtal::ImageContainerBySTLVector< DGtal::Z2i::Domain, unsigned char> Image2D;
  typedef Image2D::Domain Domain;

  Image2D image = DGtal::GenericReader<Image2D>::import(imageFilename); 
     
  
  // 2. Define a domain, a digital set, and build a digital set from an image and a predicate (threshold values)
   Image2D::Domain wholeDomain(image.domain().lowerBound()+Z2i::Point(-2,-2), image.domain().upperBound()+Z2i::Point(2,2)); 
   
   /*************************************************************/
   // Optional: if the image is too big, use the following subsampling to work on a smaller one
   // Subsample to twice smaller grid.

   //#ifdef SUBSAMP
   std::cout << "Subsample " << subsamp << std::endl;
   
   //std::vector<Image2D::Domain::Size> aGridSize2D(SUBSAMP,SUBSAMP);
   std::vector<Image2D::Domain::Size> aGridSize2D;
   aGridSize2D.push_back(subsamp);
   aGridSize2D.push_back(subsamp);
   DGtal::functors::BasicDomainSubSampler<Image2D::Domain> subSampler2D(wholeDomain, aGridSize2D, Z2i::Point(0 ,0));
   
   typedef ConstImageAdapter<Image2D,  Image2D::Domain, 
			     functors::BasicDomainSubSampler<Image2D::Domain>,  
			     Image2D::Value,
			     functors::Identity > SubsampledImage;
   
   functors::Identity df;
   // Get the new domain produces from the subsampler and define the ConstImageAdapter:
   Image2D::Domain subSampledDomain2D  = subSampler2D.getSubSampledDomain();
   SubsampledImage subsampledImage2D (image, subSampledDomain2D, subSampler2D, df);
   
   //GenericWriter<SubsampledImage>::exportFile("subsampled.pgm", subsampledImage2D );

   Image2D::Domain domain(subSampledDomain2D.lowerBound()+Z2i::Point(-2,-2),subSampledDomain2D.upperBound()+Z2i::Point(2,2));
   
 
   // 2. Define a domain, a digital set, and build a digital set from an image and a predicate (threshold values)

   typedef Z2i::DigitalSet DigitalSet;
   DigitalSet mySet(domain);
   typedef Z2i::Point Point;
   
   SetFromImage<DigitalSet>::append<SubsampledImage>(mySet, subsampledImage2D, 0 ,255);
   
   // 3. EPS output of the digital set (using Board2D)
   Board2D board;
      
   board << domain;
   board << mySet;
   board.saveEPS("digitalSet.eps");
   

   // 1. Choose an adjacency, define a Digital Object and initialize it with the digital set;

   typedef Z2i::Object4_8 DigitalObject;
   DigitalObject myObject(Z2i::dt4_8,mySet);

   vector<DigitalObject> objects;
   back_insert_iterator< vector<DigitalObject> > it( objects );
   myObject.writeComponents( it ); 

   int maxInd = 0;
   int j = 0;
   for(auto const & o : objects)
     {
       if(o.size() > objects[maxInd].size())
	 maxInd = j;
       j++;
     }

   DigitalObject obj(objects[maxInd]);
   
   
   /*************** Compute set of enclosing balls ***********************/
   // Compute the set of points not in the digital set and $4$ connected to it.
   
   DigitalSet complementSet(domain);
   
   for(auto const &p : domain)
     {
       if(!((obj.pointSet())(p)))
	 complementSet.insert(p);
     }

   typename DigitalObject::ComplementObject complement(Z2i::dt8_4,complementSet);
   typename DigitalObject::ComplementObject border = complement.border();


   board.clear();
   board << domain;
   board << border;
   board.saveEPS("outBorder.eps");
   
   
   // Compute the voronoi diagram of the out border --> use CGAL

   typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt Kernel;
   typedef CGAL::Delaunay_triangulation_2<Kernel>                                    DT;
   typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
   typedef CGAL::Identity_policy_2<DT,AT> AP; // adaptation policy --> here keep everything
   typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VoronoiD;

   typedef CGAL::Polygon_2<Kernel> Polygon_2;
   typedef typename Kernel::Point_2 Point_2;
   typedef CGAL::Circle_2<Kernel> Circle_2;
   
   VoronoiD voronoi;
   
   for(auto const &p : border)
     {
       typename VoronoiD::Site_2 s = Kernel::Point_2(p[0],p[1]);
       voronoi.insert(s);
     }

   board.clear();
   board << domain;
   board << obj;
   board.setPenColor( DGtal::Color::Red );
   int i = 0;
   for(VoronoiD::Edge_iterator ei = voronoi.edges_begin() ; ei != voronoi.edges_end() ; ei++)
     {
       if(ei->is_segment())
	 {
	   i++;
	   typename VoronoiD::Point_2 p = ei->source()->point(), q = ei->target()->point();  
	   board.drawLine(p.x().floatValue(),p.y().floatValue(),q.x().floatValue(),q.y().floatValue());
	 }
     }
   board.saveEPS("voronoi.eps");
   std::cout << "number of voronoi segments = " << i << std::endl;
   
   
   // Compute boundary of the set --> used to crop the Voronoi diagram

   std::cout << std::endl ; trace.beginBlock("Compute polygon used to crop VD");
   Z2i::KSpace ks;
   ks.init(domain.lowerBound(), domain.upperBound(),false);
   
   SurfelAdjacency<2> sAdj(true);
   
   std::vector< std::vector<Z2i::KSpace::SCell> > SCellContours;
   std::vector<std::vector<Z2i::Point>> PointsContours;
   try
     {
       // Surfaces<Z2i::KSpace>::extractAll2DSCellContours(SCellContours, ks, sAdj, obj.pointSet());
       Surfaces<Z2i::KSpace>::extractAllPointContours4C(PointsContours, ks, obj.pointSet(), sAdj);
     }
   catch(DGtal::InputException i)
     {
       trace.emphase() << "could not find a starting bel" << std::endl;
     }
   
      

  // Compute polygon used to crop the VD
  Polygon_2 P;
  
  std::vector<Z2i::Point> contour = PointsContours[0];
  
  for(auto const & p : contour)
    P.push_back(Point_2(p[0]-0.5,p[1]+0.5));

  trace.endBlock();
  
  board.clear();
  board << domain;
  board << obj;
  board.setPenColor( DGtal::Color::Red );
  for(Polygon_2::Edge_const_iterator ei = P.edges_begin(); ei != P.edges_end() ; ei++)
    {
      typename Polygon_2::Segment_2 seg = *ei;
      //std::cout << seg <<  std::endl;
      board.drawLine(seg[0].x().floatValue(),seg[0].y().floatValue(),seg[1].x().floatValue(),seg[1].y().floatValue());
    }

  board.saveEPS("polygon.eps");
  
  //*************** Compute cropped Voronoi diagram and its graph *****************//

  std::cout << std::endl ; trace.beginBlock("Compute cropped VD + graph");

  // Define graph type
  struct VertexProperty 
  { 
    Point_2 id;
    Point_2 site; // site closest to the vertex geometric embedding --> enables to retrieve the maximal ball radius. 
  };


  // Define graph type
  struct EdgeProperty 
  {
    std::set<typename Kernel::FT> listCriticalLambdas;
  };
  
   
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty> Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex_t;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge_t;
  
  CroppedVoronoiDiagram<Kernel,Graph> CroppedVD(voronoi, P);
  
  std::vector<typename Polygon_2::Segment_2> croppedVD = CroppedVD.getCroppedVoronoiDiagram();
  CroppedVD.computeGraph();
  
  trace.endBlock();

  board.clear();
  board << domain;
  board << obj;
  board.setPenColor( DGtal::Color::Red );
  for(auto const & seg : croppedVD)
    board.drawLine(seg[0].x().floatValue(),seg[0].y().floatValue(),seg[1].x().floatValue(),seg[1].y().floatValue());
  
  board.saveEPS("croppedVD.eps");
  
  
  board.clear();
  board << domain ;
  board << mySet;
  CroppedVD.drawVDGraph(board);
  board.saveEPS("VDGraph.eps");

  
  
   // //************* Compute topological ordering on VD vertices ********//

  std::cout << std::endl ; trace.beginBlock("Compute topological ordering");
  
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS, VertexProperty, EdgeProperty> DirectedGraph;
   
   std::vector<typename DirectedGraph::vertex_descriptor> topologicalOrdering;
   
   std::vector<typename DirectedGraph::edge_descriptor> topologicalOrderingEdges;

   Graph g = CroppedVD.getVDGraph();

   DirectedGraph dag = directedGraphFromTree<Graph,DirectedGraph>(g, topologicalOrdering,topologicalOrderingEdges);

   trace.endBlock();

      
   
   //*************** Compute critical balls on critical pencils ********************//
   
   std::cout << std::endl ; trace.beginBlock("Compute ordered set of critical balls on critical pencils");

   std::vector<Circle_2> orderedCriticalBalls;

   
   typedef ImageContainerBySTLVector<Domain, Circle_2> CriticalBalls;
   CriticalBalls criticalBalls(domain);
   
   std::set<typename Kernel::FT> listCriticalLambdas;

   int n = topologicalOrdering.size()-1;

   i=0;
   
   for(std::vector<typename DirectedGraph::edge_descriptor>::reverse_iterator it = topologicalOrderingEdges.rbegin() ; it != topologicalOrderingEdges.rend() ; it++)
     {
       std::cout << i << "/" << n << "\r";
       std::cout.flush();
       i++;
       

       Circle_2 b1,b2,b;
       pencilFromEdge<Kernel,DirectedGraph>(*it,dag,b1,b2);
       
       // compute the set of digital points that are in the interior of b1 but
       // not in the interior of b2. The points for which a critical
       // ball may be found on the pencil [b1,b2] are among these points 
       DigitalSet digitBall1 = digitize<Kernel,DigitalSet>(b1);
       
       DigitalSet digitBall2 = digitize<Kernel,DigitalSet>(b2);
       typename Kernel::FT lambda;
       for(auto const &p : digitBall1)
	 {
	   if(!digitBall2(p))
	     {
	       Point_2 pp(p[0],p[1]);
	       lambda = computeLambdaOfCriticalBall<Kernel,DigitalSet>(pp,b1,b2);
	       
	       if(lambda < 0 || lambda > 1)
		 std::cerr <<  "Lambda out of [0,1] interval !\n";
	       
	       b = ballOnPencil<Kernel>(lambda,b1,b2); 		   
	       
	       criticalBalls.setValue(p,b);
	       
	       dag[*it].listCriticalLambdas.insert(lambda);
	     }
	 }
     }
   
   // for all the points in the root ball, set criticalBall to the right value = root ball
   std::vector<typename DirectedGraph::edge_descriptor>::iterator itroot = topologicalOrderingEdges.begin();
   
   typename DirectedGraph::vertex_descriptor vroot = boost::target(*itroot,dag);
   Circle_2 b = circleFromVertex<Kernel>(vroot,dag);
   DigitalSet digit = digitize<Kernel,DigitalSet>(b);
   for(auto const & p : digit)
     criticalBalls.setValue(p,b);
   
   trace.endBlock();

   //*************** Compute optimal cover ********************//
   
   
   trace.beginBlock("Compute optimal cover");
   
   DigitalSet unmarked(obj.pointSet());
   std::vector<Circle_2> optimalCover;

   i=0;
   
   for(std::vector<typename DirectedGraph::edge_descriptor>::reverse_iterator it = topologicalOrderingEdges.rbegin() ; it != topologicalOrderingEdges.rend() ; it++)
     {
       std::cout << i << "/" << n << "\r";
       std::cout.flush();
       i++;
       
       Circle_2 b1,b2,b;
       pencilFromEdge<Kernel,DirectedGraph>(*it,dag,b1,b2);

       for(std::set<Kernel::FT>::iterator itL = dag[*it].listCriticalLambdas.begin() ; itL  != dag[*it].listCriticalLambdas.end() ; itL++)
	 {
	   //std::cout << *itL << "  ";
	   Circle_2 b = ballOnPencil<Kernel>(*itL,b1,b2);

	   DigitalSet digit = digitize<Kernel,DigitalSet>(b);
	   bool found = false;
	   
	   for(auto const & p : digit)
	     {
	       if(unmarked(p) && b == criticalBalls(p))
		 {
		   found = true;
		   optimalCover.push_back(b);
		   break;
		 }
	     }
	   if(found) // a ball has been added, update unmarked set
	     for(auto const & p : digit)
	       if(unmarked(p))
		 unmarked.erase(p);
	 }
     }

   // if unmarked points remain, add the root ball
   if(unmarked.size() !=0)
     {
       std::cout << "unmarked vertices remain\n";
       std::vector<typename DirectedGraph::edge_descriptor>::iterator it = topologicalOrderingEdges.begin();

       typename DirectedGraph::vertex_descriptor vroot = boost::target(*it,dag);
       Circle_2 b = circleFromVertex<Kernel>(vroot,dag);
       optimalCover.push_back(b);
       DigitalSet digit = digitize<Kernel,DigitalSet>(b);
       for(auto const & p : digit)
	 if(unmarked(p))
	   unmarked.erase(p);       
     }
   
   
   std::cout << "number unmarked = " << unmarked.size() << std::endl;
   std::cout << "Number of balls = " << optimalCover.size() << std::endl;

   trace.endBlock();


   //*************** Output results ********************//
   
   board.clear();
   board << domain;

   for(auto const & p :obj)
     board << CustomStyle( p.className(), 
			   new CustomPen( Color(180,180,180), Color(210,210,210), 1.0)) << p;
   
   board.setLineWidth(2.0);
   board.setPenColor( DGtal::Color::Green );
   board.setPenColorRGBi(50,205,50 );
   for(auto const & seg : croppedVD)
     board.drawLine(seg[0].x().floatValue(),seg[0].y().floatValue(),seg[1].x().floatValue(),seg[1].y().floatValue());
   board.setPenColor(DGtal::Color::Red);
   
   board.setFillColor(DGtal::Color::None);
   for(auto const & b : optimalCover)
     board.drawCircle(b.center().x().floatValue(), b.center().y().floatValue(), sqrt(b.squared_radius().floatValue()));
  

   board.saveEPS("optimalCover.eps");
		      
		           
}

