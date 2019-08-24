/*
 * Prim.cpp
 *
 *  Created on: 21 Aug 2019
 *      Author: Ian MacDonald
 *      
 *      
 */

#include <chrono>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <string>

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
    bool add_neighbour(int vertex, const EdgeDataType &data=EdgeDataType()) {
      neighbours.push_back(Neighbour<EdgeDataType>(vertex,data));
      return true;
    }
    bool add_neighbour(int vertex, std::istream & in) {
      EdgeDataType data;
      const bool z=in>>data;
      if (z)
        neighbours.push_back(Neighbour<EdgeDataType>(vertex,data));
      return z;
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
    bool add_neighbour(int vertex, const EdgeDataType &data=EdgeDataType()) {
      neighbours[vertex]=Neighbour<EdgeDataType>(vertex,data);
      return true;
    }
    bool add_neighbour(int vertex, std::istream & in) {
      EdgeDataType data;
      const bool z=in>>data;
      if (z)
        neighbours[vertex]=Neighbour<EdgeDataType>(vertex,data);
      return z;
    }
    // Find neighbour and return pointer to data
    const EdgeDataType *find_neighbour(int i) {
      if (neighbours[i].get_vertex()>=0) return &neighbours[i].get_data();
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
      return const_iterator(it,neighbours.cend());
    };
    const_iterator end() const {return const_iterator(neighbours.cend(),neighbours.cend());};
    const_iterator cbegin() const {
      auto it=neighbours.cbegin();
      // Find first entry
      while ((it!=neighbours.cend())&&(it->get_vertex()<=0)) ++it;
      return const_iterator(it,neighbours.cend());
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
    Graph(int nvert) : nvert(nvert), g(nvert, NeighSetType(nvert)) {};

    // Constructor to read in graph from given file name
    Graph(std::string filename) : nvert(0), g(0, NeighSetType(0)){
      std::ifstream in(filename);
      in>>nvert;
      g.resize(nvert,NeighSetType(nvert));
      int i,j;
      do {
        in>>i>>j;                 // Read edge vertices
        g[i].add_neighbour(j,in); // Add edge and get any edge data from stream in
      } while( in );
      //} while( in  && std::cout<<"> "<<i<<" "<<j<<" "<<*(g[i].find_neighbour(j))<<"\n");
      if (!in.eof()) {
        std::cout<<"Error reading file: "<<filename<<"\n";
        exit (EXIT_FAILURE);
      }
      in.close();
    }

    NeighSetType &operator[](int i){return g[i];};

    const NeighSetType &operator[](int i) const {return g[i];};

    int get_nvert() const {return nvert;};

    friend std::ostream &operator<< <>(std::ostream &out, const Graph<NeighSetType> &g);

    // Print adjacency matrix
    void print_adjacency_matrix(std::ostream &out) {
      for (int i=0; i<nvert; i++) {
        out<<i<<": ";
        for (int j=0; j<nvert; j++) {
          auto ed=g[i].find_neighbour(j);
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

template<class GraphType>
class RandomUndirectedGraph : public GraphType {
  public:
    template<class RandEdgeFunc>
    RandomUndirectedGraph(int nvert,
        std::default_random_engine &gen, // Reference to a default random engine
        double p, // probabily that and vertex-vertex connection will exist
        RandEdgeFunc ran_edge_data) // Function to create random edge data
                                    // i.e. EdgeDataType(EdgeDataType &)
        : GraphType(nvert) {

      // Create uniform distribution [0,1]
      auto dist=std::uniform_real_distribution<double>(0.0,1.0);
      // Loop over possible vertex pairs
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
// Computes a Minimum Spanning Tree (MST) for the Graph<NeighSetType>
// using Prim's algorithm.
// Returns the number of vertices in the closed set. If this is not equal
// to the total number of vertices, then the graph is not fully connected
// and the algorithm has FAILED.
// On SUCCESS:
// A MST consists of the edges (i,parent[i]) i=1,....,nvert-1
// length[i] is the length of edge (i,parent[i])
// mst_len=SUM(length[i], i=1,....,nvert-1) is the total length of the MST
//***************************************************************************
template<class NeighSetType, class DistType, class EdgeLengthFunc>
int Prim_min_span_tree(
    const graph::Graph<NeighSetType> &g,
    EdgeLengthFunc edge_data_to_length, // Function to map edge data to edge length
                                        // i.e. DistType(const EdgeDataType &)
    DistType len_large,                 // Length larger than any possible path length.
    DistType &mst_len,                  // RETURNS the length of minimum spanning tree, provided
                                        // that the graph is fully connected
    // The first nvertex entries of the following vectors are set by this routine.
    std::vector<int>      &parent,      // parent[i] RETURNS the parent vertex of vertex i.
    std::vector<DistType> &distance,    // distance[i] RETURNS the length of edge (i,parent[i]).
    std::vector<bool>     &zclosed      // zclosed[i] RETURNS whether vertex i is in the closed set.
    ) {
    int nvert=g.get_nvert();
    // Initialise the empty closed set
    int closed_size=0;
    // Until a vertex i is in the closed set, distance[i] will store the edge length
    // to the nearest neighbour in the closed set, or "len_large" if no neighbours are in
    // the closed set. parent[i] will store the nearest closed neighbour (-1 if none).
    // Once a new vertex has been added to the closed set, must update distance and
    // parent values for its non-closed neighbours.
    for (int i=0; i<nvert; i++) {
      zclosed[i]=false;
      distance[i]=len_large;
      parent[i]=-1;
    }
    int cvert; // Current vertex
    DistType clen;
    mst_len=0.0;
    using pair=std::tuple<DistType,int>;
    // Create a priority queue of pairs, such that lowest distance
    // is always at the top
    std::priority_queue<pair,std::vector<pair>,std::greater<pair>> pq;
    // Add vertex 0 with zero distance to priority queue, so that this becomes
    // the first member of the closed set.
    pq.push(std::make_pair(0.0,0));
    distance[0]=0.0;
    while(!pq.empty()) { // Loop until priority queue is empty
      std::tie(clen,cvert)=pq.top();  // Get lowest distance (top) from priority queue.
      pq.pop();                       // Remove the top pair from queue.
      if (!zclosed[cvert]) {          // Skip if already in closed set
        // Add vertex to close set
        zclosed[cvert]=true;
        closed_size++;
        mst_len+=clen;
        if (closed_size==nvert) break; // Finish if all vertices are in the closed set.

        // Update distance to closed set for non-closed neighbours
        for (auto e: g[cvert]) { // Loop over neighbours of current vertex
          const int neigh=e.get_vertex();
          if (!zclosed[neigh]) {
            const DistType neigh_dist=edge_data_to_length(e.get_data()); // Get edge length
            if (neigh_dist>=0.0 && neigh_dist<distance[neigh]) {
              // Distance is less than previous distance, add this distance to priority queue.
              // Note: this vertex may already have greater distances in the priority queue.
              // However, these will be ignored once the vertex is in the closed set.
              pq.push(std::make_pair(neigh_dist,neigh));
              // Set parent vertex and distance to parent
              parent[neigh]=cvert;
              distance[neigh]=neigh_dist;
            }
          }
        }
      }
    }
    if (closed_size<nvert) mst_len=len_large; // Check graph is fully connected.
    return closed_size;
}

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

//***************************************************************************
// Main routine reads creates a graph using a
// The Prim_min_span_tree routine is called for this graph and the
// Minimum spanning tree (MST) information is written out.
//***************************************************************************

enum class Colour { RED, GREEN, BLUE, First=RED, Last=BLUE };
const std::vector<std::string> cmap{"RED","GREEN","BLUE"};
inline std::ostream &operator<<(std::ostream &out,
    const Colour colour) {
    out<<cmap[static_cast<int>(colour)];
    return out;
}

struct edge_data {
    int weight;
    Colour colour;
    edge_data(int weight, Colour colour=Colour::RED): weight(weight), colour(colour) {};
};
inline std::ostream &operator<<(std::ostream &out,
    edge_data const data) {
    out<<data.weight<<" "<<data.colour;
    return out;
}

int main() {
  const bool zverbose=false;

  // Use int distances
  using DistType=int;
  using EdgeDataType=edge_data; // The edge data consists of the same type

  // Create graph from given file name.
  //graph::Graph<graph::NeighbourSetAdjList<EdgeDataType>> gr("SampleTestData_mst_data.dat");
  //graph::Graph<graph::NeighbourSetAdjMat<EdgeDataType>> gr("SampleTestData_mst_data.dat");

  // Create a random Graph

  // Get random seed
  //unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed=12345678; // Fixed seed for run-to-run comparisons
  std::cout<<"seed="<<seed<<"\n";
  // Create random generator initialised from seed
  std::default_random_engine gen(seed);

  const DistType wmin=1, wmax=10;
  const double prob=0.8;
  const int nvert=20;

  using GraphType=graph::Graph<graph::NeighbourSetAdjList<EdgeDataType>>;
  graph::RandomUndirectedGraph<GraphType> gr(
      nvert, gen, prob,
      // Lambda to compute a random edge distance in the range w1<=w<=w2
      [wmin,wmax](std::default_random_engine &gen)
      {
    auto dist1=get_right_distribution<DistType>(wmin,wmax);
    auto dist2=get_right_distribution<int>(static_cast<int>(Colour::First),
                                           static_cast<int>(Colour::Last));
    return edge_data(dist1(gen),static_cast<Colour>(dist2(gen)));
      }
  );

  // Write out graph
  std::cout<<"nvert="<<gr.get_nvert()<<"\n";
  // Ouput edges and data for graph
  std::cout<<"Adjacency List:\n"<<gr<<"\n";
  if (zverbose) {
    // Output adjacency matrix for graph
    std::cout<<"Adjacency Matrix:\n";
    gr.print_adjacency_matrix(std::cout);
  }

  const DistType len_large=1000000000; // Larger than any possible path length
  DistType mst_len;
  std::vector<int>      parent(nvert); // (edges (i,parent[i]) i=1,i...,nvert-1
                                       // will be the minimum spanning tree
  std::vector<DistType> length(nvert); // length[i] will be the length of edge (i,parent[i])
  std::vector<bool>     zclosed(nvert);// Stores closed set membership

  for (int ic=static_cast<int>(Colour::First); ic<=static_cast<int>(Colour::Last); ic++) {
    const Colour colour=static_cast<Colour>(ic);
    std::cout<<"Colour="<<colour<<"\n";

    // Compute MST using Prim's algorithm`
    const int nclosed=Prim_min_span_tree(gr,
        [colour](const EdgeDataType &data){
      return static_cast<DistType>(data.colour==colour ? data.weight : -1);
    },
    len_large,
    mst_len,parent,length,zclosed);

    if (zverbose) {
      // Print the returned data vectors
      std::cout<<"nclosed="<<nclosed<<"\n";
      std::cout<<"         i: ";
      for (int i=1; i<nvert; i++) std::cout<<std::setw(5)<<i<<" ";
      std::cout<<"\n";
      std::cout<<" parent[i]: ";
      for (int i=1; i<nvert; i++) std::cout<<std::setw(5)<<parent[i]<<" ";
      std::cout<<"\n";
      std::cout<<" length[i]: ";
      for (int i=1; i<nvert; i++) std::cout<<std::setw(5)<<length[i]<<" ";
      std::cout<<"\n";
      std::cout<<"zclosed[i]: ";
      for (int i=1; i<nvert; i++) std::cout<<std::setw(5)<<zclosed[i]<<" ";
      std::cout<<"\n";
    }

    if (nclosed==nvert) {
      // Output only valid is all vertices are in the closed set

      std::cout<<"Minimum spanning tree length= "<<mst_len<<"\n";

      // Output edges of MST
      int count=0;
      DistType len=0; // Check length of MST
      for (int i=1; i<nvert; i++) {
        //std::cout<<" Edge "<<parent[i]<<"-"<<i<<" Weight "<<length[i]<<"\n";
        std::cout<<" Edge "<<parent[i]<<"-"<<i<<" Weight "<<length[i]<<": "
            <<*(gr[i].find_neighbour(parent[i]))<<"\n";
        count++;
        len+=length[i];
      }
      std::cout<<"count="<<count<<" length="<<len<<"\n";
    } else {
      std::cout<<"Graph not full connected\n";
    }

  }

  return 0;
}
