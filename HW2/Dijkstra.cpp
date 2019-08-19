/*
Code Output:
average shortest distance for 50 vertices with probability 0.2 and weight in (1,10) is 6.97778 (99898 samples)
average shortest distance for 50 vertices with probability 0.4 and weight in (1,10) is 4.7069 (100000 samples)


There are a number approaches for implementing a graph ADT. I wanted to develop a framework that would be efficient for the Dijkstra method, but could also be tailored to be efficient for other algorithms. The key requirement for Dijkstra is to loop over all vertices connected to a particular vertex, obtaining the weights for each of the corresponding edges. I will refer to the connected vertices as NEIGHBOURS.

My framework is based on a neighbour set container class that holds the set of neighbours (with templated edge data) for a single vertex. An iterator is provided to make it easy to loop over the neighbours. My "Graph" class is now just a container for the neighbour sets of all the vertices and is templated on the neighbour set container class.

The interface for the neighbour set container is independent of its internal representation. In regards to Dijkstra, it seems logical that the Adjacency List Scheme is the most efficient. The class "NeighbourSetAdjList" stores a std::vector of the neighbours of a vertex. This approach makes finding whether a particular vertex is a neighbour ("find_neighbour" method) inefficient, but this is not required by Dijkstra.

I have also implemented the class "NeighbourSetAdjMat" that uses the alternative Adjacency Matrix strategy. Here a std::vector stores a potential neighbour entry for all vertices. This makes the find operation efficient, but looping over neighbours less efficient than the Adjacency List strategy. I would expect the Adjacency Matrix approach to be more costly than the Adjacency List one, because of the greater inefficiency in looping over vertex neighbours for Dijkstra. As the edge probability increases, I would expect the cost of the Adjacency Matrix strategy to approach that of Adjacency List, but still be more costly. However, the actual results are slightly counter-intuitive, in that the Adjacency Matrix formulation actually becomes cheaper than the Adjacency List one for large enough edge probabilities. For example, this is at a probability of only about 0.2 with 50 vertices. The switch over probability increases with increasing number of vertices. It seems likely that cache efficiency may explain why the Adjacency Matrix approach becomes cheaper. It has the advantage that all the vectors are of known size when the "Graph" constructor is called. This would allow them to be allocated consecutively in memory.

The "Graph" class is extended to the "Dijkstra" class in order to supply the method to carry out the Dijkstra shortest method algorithm. The "Graph" class, or any class derived from it (e.g. "Dijkstra") can be extended to the "RandomUndirectedGraph" class. Its constructor sets up a undirected random graph, with a specified edge probability and random edge weights, in the range wmin<=w<=wmax.

The main program asks the user for the model parameters. It then computes the average shortest path length over the requested number of randomly generated graphs. This is done for both the Adjacency List and Adjacency Matrix based neighbour set containers, allowing the computational costs to be compared.

Observations from results:
For a fixed number of vertices, the mean shortest path reduces with increasing edge probability. This is because there is an increasing probability of finding a shorter path between vertices.

For a fixed probability, the mean shortest path reduces with increasing number of vertices

C++ techniques learned:
1) How to provide iterators for a custom container class
2) How to build classes using templates
3) Creating and passing lambdas to methods
4) How to use the new STL random number generation
*/

/*
 * Dijkstra.cpp
 *
 *  Created on: 13 Aug 2019
 *      Author: Ian MacDonald
 *      
 *      
 */

#include <chrono>
#include <cstddef>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace graph {

//***************************************************************************
// A graph is described by the neighbours connected to each vertex.
// The class "Neighbour<EdgeDataType>" is used to describe a
// single neighbour and the data associated with the corresponding edge.
//***************************************************************************
template<class EdgeDataType> class Neighbour;
template<class EdgeDataType> std::ostream &operator<<(std::ostream &out,
    Neighbour<EdgeDataType> const neigh);

template<class EdgeDataType>
class Neighbour {
  public:
    // Constructor
    Neighbour() : vertex(-1), data() {}
    Neighbour(int vertex, EdgeDataType const &data): vertex(vertex), data(data) {};
    // Get vertex number
    int get_vertex() const {return vertex;};
    // Get reference to edge data
    EdgeDataType &get_data() {return data;};
    // Get const reference to edge data
    const EdgeDataType &get_data() const {return data;};

    friend std::ostream &operator<< <>(std::ostream &out,
        Neighbour<EdgeDataType> const neigh);
  private:
    int vertex; // Connecting vertex
    EdgeDataType data;     // Edge data
};

template<class EdgeDataType> inline std::ostream &operator<<(std::ostream &out,
    Neighbour<EdgeDataType> const neigh) {
    out<<"("<<neigh.vertex<<" "<<neigh.data<<") ";
    return out;
}

//***************************************************************************
// The class "NeighbourSetAdjList<EdgeDataType>" holds the set of neighbours 
// for a particular vertex using the Adjacency List approach.
// This implementation stores a std::vector of the "Neighbour<EdgeDataType>" objects.
// Iterators are provided to allow the the user to iterate over the set of neighbours
//***************************************************************************
template<class EdgeDataType> class NeighbourSetAdjList;
template<class EdgeDataType> class NeighbourSetAdjList_iterator;
template<class EdgeDataType> class NeighbourSetAdjList_const_iterator;
template<class EdgeDataType> std::ostream &operator<<(std::ostream &out,
    NeighbourSetAdjList<EdgeDataType> const list);

template<class EdgeDataType>
class NeighbourSetAdjList {
  public:
    using iterator=NeighbourSetAdjList_iterator<EdgeDataType>;
    using const_iterator=NeighbourSetAdjList_const_iterator<EdgeDataType>;

    // Constructor
    NeighbourSetAdjList(int nvert){}; //nvert not needed for this implementation

    // Add neighbour
    void add_neighbour(int vertex, const EdgeDataType &data) {
      neighbours.push_back(Neighbour<EdgeDataType>(vertex,data));
    }

    // Find neighbour and return pointer to data (nullptr if not found)
    const EdgeDataType *find_neighbour(int i) {
      for (auto &neigh: neighbours) {
        if (i==neigh.get_vertex()) return &neigh.get_data();
      }
      return nullptr; // Neighbour i not found
    }
    // iterator begin and end
    iterator begin() {return iterator(neighbours.begin());};
    iterator end() {return iterator(neighbours.end());};
    const_iterator begin() const {return const_iterator(neighbours.cbegin());};
    const_iterator end() const {return const_iterator(neighbours.cend());};
    const_iterator cbegin() const {return const_iterator(neighbours.cbegin());};
    const_iterator cend() const {return const_iterator(neighbours.cend());};

    friend std::ostream &operator<< <>(std::ostream &out,
        NeighbourSetAdjList<EdgeDataType> const list);

  private:
    // Vector of neighbours
    std::vector<Neighbour<EdgeDataType>> neighbours;
};

template<class EdgeDataType>
inline std::ostream &operator<<(std::ostream &out,
    NeighbourSetAdjList<EdgeDataType> const list) {
    for (auto neigh: list.neighbours)
      out<<neigh<<" ";
    return out;
}

//*************************************************
// iterator class for "NeighbourSetAdjList"
//*************************************************
// iterator class
template<class EdgeDataType>
class NeighbourSetAdjList_iterator {
  public:
    using self_type=NeighbourSetAdjList_iterator;
    using iterator_category=std::forward_iterator_tag;
    using value_type=Neighbour<EdgeDataType>;
    using difference_type=std::ptrdiff_t;
    using pointer=Neighbour<EdgeDataType>*;
    using reference=Neighbour<EdgeDataType>&;
    // Constructor
    NeighbourSetAdjList_iterator(
        typename std::vector<value_type>::iterator const &it): it(it){};
    reference operator*() const {return *it;};
    pointer operator->() const {return it;};
    self_type& operator++() {it++; return *this;};
    self_type operator++(int) {
      self_type tmp(*this);
      it++;
      return tmp;
    }
    bool operator==(const self_type &rhs) const {return it==rhs.it;};
    bool operator!=(const self_type &rhs) const {return it!=rhs.it;};
  private:
    typename std::vector<value_type>::iterator it;
};

//*************************************************
// const_iterator class for "NeighbourSetAdjList"
//*************************************************
template<class EdgeDataType>
class NeighbourSetAdjList_const_iterator {
  public:
    using self_type=NeighbourSetAdjList_const_iterator;
    using iterator_category=std::forward_iterator_tag;
    using value_type=Neighbour<EdgeDataType>;
    using difference_type=std::ptrdiff_t;
    using pointer=const Neighbour<EdgeDataType>*;
    using reference=const Neighbour<EdgeDataType>&;
    // Constructor
    NeighbourSetAdjList_const_iterator(
        typename std::vector<value_type>::const_iterator const &it): it(it){};
    reference operator*() const {return *it;};
    pointer operator->() const {return it;};
    self_type& operator++() {it++; return *this;};
    self_type operator++(int) {
      self_type tmp(*this);
      it++;
      return tmp;
    }
    bool operator==(const self_type &rhs) const {return it==rhs.it;};
    bool operator!=(const self_type &rhs) const {return it!=rhs.it;};
  private:
    typename std::vector<value_type>::const_iterator it;
};

//***************************************************************************
// The class "NeighbourSetAdjMat<EdgeDataType>" holds the set of neighbours 
// for a particular vertex using the Adjacency Matrix approach.
// This implementation stores a std::vector with an entry for the all possible neighbours.
// Only the entries corresponding to actual neighbours have valid vertex numbers set.
// Iterators are provides to allow the the user to iterate over the set of neighbours
// The iterators are more complicated and inefficient for this implementation,
// because they need to skip over the entries that are not neighbours.
//***************************************************************************
template<class EdgeDataType> class NeighbourSetAdjMat;
template<class EdgeDataType> class NeighbourSetAdjMat_iterator;
template<class EdgeDataType> class NeighbourSetAdjMat_const_iterator;
template<class EdgeDataType> std::ostream &operator<<(std::ostream &out,
    NeighbourSetAdjMat<EdgeDataType> const list);

template<class EdgeDataType>
class NeighbourSetAdjMat {
  public:
    using iterator=NeighbourSetAdjMat_iterator<EdgeDataType>;
    using const_iterator=NeighbourSetAdjMat_const_iterator<EdgeDataType>;
    // Constructor
    NeighbourSetAdjMat(int nvert) : neighbours(nvert) {};
    void add_neighbour(int vertex, const EdgeDataType &data) {
      neighbours[vertex]=Neighbour<EdgeDataType>(vertex,data);
    }
    // Find neighbour and return pointer to data
    const EdgeDataType *find_neighbour(int i) {
      auto neigh=neighbours[i];
      if (neigh.get_vertex()>0) return &neigh.get_data();
      return nullptr; // Neighbour i not found
    }
    // iterator begin and end
    iterator begin() {
      auto it=neighbours.begin();
      // Find first entry
      while ((it!=neighbours.end())&&(it->get_vertex()<=0)) ++it;
      return iterator(it,neighbours.end());
    };
    iterator end() {return iterator(neighbours.end(),neighbours.end());};
    const_iterator begin() const {
      auto it=neighbours.cbegin();
      // Find first entry
      while ((it!=neighbours.cend())&&(it->get_vertex()<=0)) ++it;
      return iterator(it,neighbours.cend());
    };
    const_iterator end() const {return const_iterator(neighbours.cend(),neighbours.cend());};
    const_iterator cbegin() const {
      auto it=neighbours.cbegin();
      // Find first entry
      while ((it!=neighbours.cend())&&(it->get_vertex()<=0)) ++it;
      return iterator(it,neighbours.cend());
    };
    const_iterator cend() const {return const_iterator(neighbours.cend(),neighbours.cend());};

    friend std::ostream &operator<< <>(std::ostream &out,
        NeighbourSetAdjMat<EdgeDataType> const list);

  private:
    // Vector of neighbours
    std::vector<Neighbour<EdgeDataType>> neighbours;
};

template<class EdgeDataType>
inline std::ostream &operator<<(std::ostream &out,
    NeighbourSetAdjMat<EdgeDataType> const set) {
    for (auto neigh: set.neighbours)
      if (neigh.get_vertex()>0) out<<neigh<<" ";
    return out;
}

//*************************************************
// iterator class for "NeighbourSetAdjMat"
//*************************************************
template<class EdgeDataType>
class NeighbourSetAdjMat_iterator {
  public:
    using self_type=NeighbourSetAdjMat_iterator;
    using iterator_category=std::forward_iterator_tag;
    using value_type=Neighbour<EdgeDataType>;
    using difference_type=std::ptrdiff_t;
    using pointer=Neighbour<EdgeDataType>*;
    using reference=Neighbour<EdgeDataType>&;
    // Constructor
    NeighbourSetAdjMat_iterator(
        typename std::vector<value_type>::iterator const &it,
        typename std::vector<value_type>::iterator const &itend)
    : it(it), itend(itend) {};
    reference operator*() const {return *it;};
    pointer operator->() const {return it;};
    self_type& operator++() {
      it++;
      // Find next valid neighbour
      while ((it!=itend)&&(it->get_vertex()<=0)) ++it;
      return *this;
    };
    self_type operator++(int) {
      self_type tmp(*this);
      // Find next valid neighbour
      while ((it!=itend)&&(it->get_vertex()<=0)) ++it;
      return tmp;
    }
    bool operator==(const self_type &rhs) const {return it==rhs.it;};
    bool operator!=(const self_type &rhs) const {return it!=rhs.it;};
  private:
    typename std::vector<value_type>::iterator it;
    typename std::vector<value_type>::iterator itend;
};

//*************************************************
// const_iterator class for "NeighbourSetAdjMat"
//*************************************************
template<class EdgeValueType>
class NeighbourSetAdjMat_const_iterator {
  public:
    using self_type=NeighbourSetAdjMat_const_iterator;
    using iterator_category=std::forward_iterator_tag;
    using value_type=Neighbour<EdgeValueType>;
    using difference_type=std::ptrdiff_t;
    using pointer=const Neighbour<EdgeValueType>*;
    using reference=const Neighbour<EdgeValueType>&;
    // Constructor
    NeighbourSetAdjMat_const_iterator(
        typename std::vector<value_type>::const_iterator const &it,
        typename std::vector<value_type>::const_iterator const &itend)
    : it(it), itend(itend) {};
    reference operator*() const {return *it;};
    pointer operator->() const {return it;};
    self_type& operator++() {
      it++;
      // Find next valid neighbour
      while ((it!=itend)&&(it->get_vertex()<=0)) ++it;
      return *this;
    }
    self_type operator++(int) {
      self_type tmp(*this);
      // Find next valid neighbour
      while ((it!=itend)&&(it->get_vertex()<=0)) ++it;
      return tmp;
    }
    bool operator==(const self_type &rhs) const {return it==rhs.it;};
    bool operator!=(const self_type &rhs) const {return it!=rhs.it;};
  private:
    typename std::vector<value_type>::const_iterator it;
    typename std::vector<value_type>::const_iterator itend;
};

//***************************************************************************
// The "Graph<NeighSetType>" class is the container of the neighbour sets for
// all vertices.
// The template parameter NeighSetType can either be 
// NeighbourSetAdjList<EdgeValueType> or  NeighbourSetAdjMat<EdgeValueType>,
// or some other to be defined neighbour set class.
// The [] operator is overloaded to return a reference to the neighbour set
// for the requested vertex.
//***************************************************************************
template<class NeighSetType> class Graph;
template<class NeighSetType> std::ostream &operator<<(std::ostream &out, const Graph<NeighSetType> &g);

template<class NeighSetType>
class Graph {
  public:
    //Constructor
    Graph(int nvert): nvert(nvert), g(nvert, NeighSetType(nvert)) {};

    NeighSetType &operator[](int i){return g[i];};

    const NeighSetType &operator[](int i) const {return g[i];};

    int get_nvert() const {return nvert;};

    friend std::ostream &operator<< <>(std::ostream &out, const Graph<NeighSetType> &g);

    void print_adjacency(std::ostream &out) {
      for (int i=0; i<nvert; i++) {
        out<<i<<": ";
        for (int j=0; j<nvert; j++) {
          auto *ed=g[i].find_neighbour(j);
          if (ed) {
            out<<"("<<j<<" "<<*ed<<")";
          } else {
            out<<"("<<j<<" false)";
          }
        }
        out<<"\n";
      }
    }
  protected:
    int nvert; // Number of vertices
    std::vector<NeighSetType> g;
};

template<class NeighSetType>
inline std::ostream &operator<<(std::ostream &out, const Graph<NeighSetType> &g){
    for (int i=0; i<g.nvert; i++)
      out<<i<<": "<<g[i]<<"\n";
    return out;
}

//***************************************************************************
// The class "Dijkstra<NeighSetType>" extends "Graph<NeighSetType>"
// adding a method to perform the Dijkstra shortest path algorithm.
//***************************************************************************
template<class NeighSetType>
class Dijkstra : public Graph<NeighSetType> {
  public:
    Dijkstra(int nvert): Graph<NeighSetType>(nvert) {};

    template<class DistType, class EdgeDistFunc>
    int shortest_distance(
        int start_vert,  // Vertex to start from.
        int end_vert,    // Stop when shortest distance has been found for this vertex.
                         // Setting this negative will mean that the method will only stop
                         // after the algorithm has visited all the vertices that it is
                         // possible to reach. The graphs is then fully connected
                         // only if the number of vertices visited (the return value) is
                         // equal to the total number of vertices.
        EdgeDistFunc edge_data_to_distance, // Function to map edge data to edge distance
                                            // i.e. DistType(const EdgeDataType &)
        DistType large, // Distance larger than any possible path length.
        // The following three vectors must be at least length nvertex:
        std::vector<bool> &zvisited,     // Set to true when a vertex has been visited.
        std::vector<DistType> &distance, // Shortest path distance to each vertex (if visited).
        std::vector<int> &parent         // Previous vertex for shortest distance path (if visited).
    ) {
        // Initialise return data
        int n_visited=0;
        for (int i=0; i<this->nvert; i++) {
          zvisited[i]=false;
          distance[i]=large;
          parent[i]=-1;
        }
        int cvert; // Current vertex visiting
        DistType cdist; // Shortest distance to current vertex
        // Type pair is the tuple (distance, vertex)
        using pair=std::tuple<DistType,int>;
        // Create std priority queue of pairs, such that lowest distance is
        // always at the top
        std::priority_queue<pair,std::vector<pair>,std::greater<pair>> pq;
        // Add starting vertex (zero distance) to priority queue
        pq.push(std::make_pair(0.0,start_vert));
        distance[start_vert]=0.0;
        while(!pq.empty()) { // Loop until priority queue is empty
          std::tie(cdist,cvert)=pq.top(); // Get lowest distance (top) from priority queue
          pq.pop(); // Remove the top pair from queue
          if (!zvisited[cvert]) { // Ignore if already visited
            zvisited[cvert]=true; // Record as visited
            n_visited++;
            // Finished if all vertices have been visited or current
            // vertex is end_vertex
            if ((n_visited==this->nvert) || (cvert==end_vert)) break;
            // Iterate over the neighbours of current vertex
            for (auto e: this->g[cvert]) {
              const int neigh=e.get_vertex();
              if (!zvisited[neigh]) { // Ignore if already visited neighbour
                // Distance to is neighbour is current distance + edge distance
                const DistType neigh_dist=cdist+edge_data_to_distance(e.get_data());
                //std::cout<<cvert<<" "<<neigh<<" "<<distance[neigh]<<
                //" "<<neigh_dist<<"\n";
                if (neigh_dist<distance[neigh]) {
                  // Distance is less than best current distance found to neighbour.
                  // Add this distance to priority queue.
                  // Note: this vertex may already have greater distances on the
                  // priority queue. However, these will be ignored because
                  // zvisited is true for the neighbour
                  pq.push(std::make_pair(neigh_dist,neigh));
                  // Set final (shortest) distance value and set parent vertex
                  distance[neigh]=neigh_dist;
                  parent[neigh]=cvert;
                }
              }
            }
          }
        }
        return n_visited;
    }
};

//***************************************************************************
// This class can be used to extend the "Graph<NeighSetType>" class, or
// any class derived from it, such as "Dijkstra<NeighSetType>",
// with the constructor creating a undirectional random graph with
// random edge data.
//***************************************************************************
template<class GraphType>
class RandomUndirectedGraph : public GraphType {
  public:
    template<class RandEdgeFunc>
    RandomUndirectedGraph(int nvert,
        std::default_random_engine &gen, // Reference to a default random engine.
        double p,                        // probabily that and vertex-vertex 
	                                 // connection will exist.
        RandEdgeFunc ran_edge_data)      // Function to create random edge data
                                         // i.e. EdgeDataType(EdgeDataType &)
        : GraphType(nvert) {

      // Create uniform distribution [0,1]
      auto dist=std::uniform_real_distribution<double>(0.0,1.0);
      // Loop over all possible vertex pairs
      for (int i=0; i<nvert; i++) {
        for (int j=i+1; j<nvert; j++) {
          // Generate random value and test against requested probability
          if (dist(gen)<p) {
            // Get random edge data
            const auto data=ran_edge_data(gen);
            // Add symmetric neighbour-neighbour relationships
            this->g[i].add_neighbour(j,data);
            this->g[j].add_neighbour(i,data);
          }
        }
      }
    }
};

} // namespace graph

//***************************************************************************
// This template function allows the code to automatically switch between
// std::uniform_real_distribution for real types and
// std::uniform_int_distribution for integral types
//***************************************************************************
template <typename T,
typename Distribution = typename std::conditional<
std::is_integral<T>::value,
std::uniform_int_distribution<T>,
std::uniform_real_distribution<T>>::type>
Distribution get_right_distribution(const T a, const T b) {
  return Distribution(a, b);
}

template<class NeighSetType, class EdgeDataType, class DistType>
void montecarlo(int nvert, double prob, DistType wmin, DistType wmax,
    int nsamples, int outfreq, unsigned seed, bool zverbose);

//***************************************************************************
// This routine routine performs some simple tests
//***************************************************************************
void simple_tests()
{
  unsigned seed=12345678; // Fixed seed for run-to-run comparisons

  montecarlo<graph::NeighbourSetAdjList<int>, int, int>
  (10, 0.5, 1, 10, 1, 0, seed, true);
  montecarlo<graph::NeighbourSetAdjList<double>, double, double>
  (10, 0.5, 1, 10, 1, 0, seed, true);
}

//***************************************************************************
// Main routine reads in the run parameters and then calls the montecarlo
// routine, first with the "NeighbourSetAdjList" scheme and then with the
// "NeighbourSetAdjMat" scheme. This allows a comparison of the performance of
// the two different graph storage schema
//***************************************************************************
int main() {

  // simple_tests();

  // Use double distances
  using DistType=double;  // This is the distance type
  //using DistType=int;   // Also works with int distances
  using EdgeDataType=DistType; // The edge data consists of the same type

  int nvert; // Number of vertices
  double prob; // Probability of each possible edge existing
  double wmin; // Minimum edge weight >0
  double wmax; // Maximum edge weigth >=wmin
  int nsamples; // Total number of random graphs to generate
  int outfreq; // Print out running average of shortest distance every outfreq samples
  // Ask for parameters
  std::cout<<"Enter number of vertices >";
  std::cin>>nvert;
  std::cout<<"Enter probability [0,1] >";
  std::cin>>prob;
  std::cout<<"Enter min egde weight wmin (>0) >";
  std::cin>>wmin;
  std::cout<<"Enter max edge weight wmax (>wmin) >";
  std::cin>>wmax;
  std::cout<<"Enter total number of samples >";
  std::cin>>nsamples;
  std::cout<<"Enter output frequency >";
  std::cin>>outfreq;
  // Echo parameters
  std::cout<<"nvert="<<nvert<<" prob="<<prob<<" wmin="<<wmin<<
      " wmax="<<wmax<<" nsamples="<<nsamples<<
      " outfreq="<<outfreq<<"\n";

  // Get random seed
  unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
  //unsigned seed=12345678; // Fixed seed for run-to-run comparisons
  std::cout<<"seed="<<seed<<"\n";

  std::cout<<"ADJACANCY LIST METHOD\n";
  montecarlo<graph::NeighbourSetAdjList<EdgeDataType>, EdgeDataType, DistType>
  (nvert, prob, wmin, wmax, nsamples, outfreq, seed, false);

  std::cout<<"ADJACANCY MATRIX METHOD\n";
  montecarlo<graph::NeighbourSetAdjMat<EdgeDataType>, EdgeDataType, DistType>
  (nvert, prob, wmin, wmax, nsamples, outfreq, seed, false);

  return 0;
}

//***************************************************************************
// This routine performs a Monte Carlo simulation consisting of "nsamples"
// random graphs. For each fully connected graph the average shortest path
// is computed and added to the running average.
// Finally the average shortest path is printed as well as the time taken.
//***************************************************************************
template<class NeighSetType, class EdgeDataType, class DistType>
void montecarlo(int nvert, double prob, DistType wmin, DistType wmax,
    int nsamples, int outfreq, unsigned seed, bool zverbose) {

    // Create random generator initialised from seed
    std::default_random_engine gen(seed);

    const DistType large=1000000000; // Larger than any possible path length
    const int start=0; // Start vertex
    const int end=-1;  // End vertex set to -1 so Dijkstra doesn't stop until visited all possible vertices
    double accumavg=0.0; // Running average of shortest distance
    int ntot=0; // Total number of successful (connected) samples

    // Vector data set by shortest_distance method
    std::vector<bool> zvisited(nvert); // On return will flag whether each vertex has been visited
    std::vector<DistType> distance(nvert); // On return will hold shortest distance for visited vertex
    std::vector<int> parent(nvert); // On return will hold the parent vertex for each visited vertex

    //Alias for the Dijkstra graph class
    using Dijkstra_graph=graph::Dijkstra<NeighSetType>;

    // Record the start time of simulation
    auto tstart = std::chrono::high_resolution_clock::now();

    // Loop over random graphs
    for (int count=0; count<nsamples; count++) {

      // Create a random Graph
      graph::RandomUndirectedGraph<Dijkstra_graph> rg(
          nvert, gen, prob,
          // Lambda to compute a random edge distance in the range w1<=w<=w2
          [wmin,wmax](std::default_random_engine &gen)
          {
            auto dist=get_right_distribution<DistType>(wmin,wmax);
            return static_cast<EdgeDataType>(dist(gen));
          }
      );

      if (zverbose) {
        // Ouput edges and data for graph
        std::cout<<rg<<"\n";
        // Output adjacency matrix for graph
        rg.print_adjacency(std::cout);
      }

      // Compute shortest distances using Dijkstra. returns the number of vertices visited
      const int nvisited=rg.shortest_distance(start,end,
          [](const EdgeDataType &data){return static_cast<DistType>(data);}, large,
          zvisited,distance,parent);

      if (zverbose) {
        // Print the returned data vectors
        std::cout<<"nvisited="<<nvisited<<"\n";
        for (int i=1; i<nvert; i++) std::cout<<zvisited[i]<<" ";
        std::cout<<"\n";
        for (int i=1; i<nvert; i++) std::cout<<distance[i]<<" ";
        std::cout<<"\n";
        for (int i=1; i<nvert; i++) std::cout<<parent[i]<<" ";
        std::cout<<"\n";
      }

      // If fully connected
      if (nvisited==nvert) {
        // Compute average shortest distance
        const double avgdist=std::accumulate(distance.cbegin()+1,
            distance.cend(), static_cast<DistType>(0))/(nvert-1.0);
        ntot++;
        // Update the running average
        accumavg=accumavg*((ntot-1.0)/ntot)+avgdist/ntot;
        // Print the running average every outfreq graphs
        if ((outfreq>0) && (ntot%outfreq==0))
          std::cout<<"ntot="<<ntot<<" avg. shortest dist="<<accumavg<<"\n";
      }
    }

    // Compute and print the time taken
    auto tfinish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = tfinish - tstart;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    // Print out the results
    std::cout<<nsamples-ntot<<" graphs out of "<<nsamples<<" were not fully connected\n";
    std::cout<<"average shortest distance for "<<nvert<<" vertices with probability "<<
        prob<<" and weight in ("<<wmin<<","<<wmax<<") is "<<
        accumavg<<" ("<<ntot<<" samples) \n";
}
