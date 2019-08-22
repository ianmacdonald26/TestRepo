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

    Graph(std::string filename) : nvert(0), g(0, NeighSetType(0)){
      std::ifstream in(filename);
      in>>nvert;
      g.resize(nvert,NeighSetType(nvert));
      int i,j;
      do {
        in>>i>>j;
        g[i].add_neighbour(j,in);
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

    void print_adjacency(std::ostream &out) {
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

//***************************************************************************
//***************************************************************************

template<class NeighSetType, class DistType, class EdgeDistFunc>
int min_span_tree(
    const Graph<NeighSetType> &g,
    EdgeDistFunc edge_data_to_distance, // Function to map edge data to edge distance
    // i.e. DistType(const EdgeDataType &)
    DistType large, // Distance larger than any possible path length.
    DistType &min_dist, // Length of minimum spanning tree, provided all vertices
    // have been visited (graph full connected)
    // The following three vectors must be at least length nvertex:
    std::vector<bool> &zvisited,     // Set to true when a vertex has been visited.
    std::vector<DistType> &distance, // Shortest path distance to each vertex (if visited).
    std::vector<int> &parent         // Previous vertex for shortest distance path (if visited).
) {
    int nvert=g.get_nvert();
    // Initialise return data
    int n_visited=0;
    for (int i=0; i<nvert; i++) {
      zvisited[i]=false;
      distance[i]=large;
      parent[i]=-1;
    }
    int cvert; // Current vertex visiting
    DistType cdist;
    min_dist=0.0;
    using pair=std::tuple<DistType,int>;
    // Create std priority queue of pairs, such that lowest distance
    // always at the top
    std::priority_queue<pair,std::vector<pair>,std::greater<pair>> pq;
    // Add vertex 0 with zero distance to priority queue
    pq.push(std::make_pair(0.0,0));
    distance[0]=0.0;
    while(!pq.empty()) { // Loop until priority queue is empty
      std::tie(cdist,cvert)=pq.top(); // Get lowest distance (top) from priority queue
      pq.pop(); // Remove the top pair from queue
      if (!zvisited[cvert]) { // Ignore if already visited
        zvisited[cvert]=true; // Record as visited
        n_visited++;
        min_dist+=cdist;
        // Finished if all vertices have been visited or current
        // vertex is end_vertex
        if (n_visited==nvert) break;
        // Iterate over the neighbours of current vertex
        for (auto e: g[cvert]) {
          const int neigh=e.get_vertex();
          if (!zvisited[neigh]) { // Ignore if already visited neighbour
            // Distance to is neighbour is current distance + edge distance
            const DistType neigh_dist=edge_data_to_distance(e.get_data());
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
    if (n_visited<nvert) min_dist=large;
    return n_visited;
}

} // namespace graph


//***************************************************************************
// Main routine reads in the run parameters and then calls the montecarlo
// routine, first with the "NeighbourSetAdjList" scheme and then with the
// "NeighbourSetAdjMat" scheme. This allows a comparison of the performance of
// the two different graph storage schema
//***************************************************************************
int main() {
  const bool zverbose=true;

  // Use int distances
  using DistType=int;
  using EdgeDataType=DistType; // The edge data consists of the same type

  graph::Graph<graph::NeighbourSetAdjList<EdgeDataType>> gr("SampleTestData_mst_data.dat");
  //graph::Graph<graph::NeighbourSetAdjMat<EdgeDataType>> gr("SampleTestData_mst_data.dat");
  const int nvert=gr.get_nvert();

  if (zverbose) {
    std::cout<<"nvert="<<nvert<<"\n";
    // Ouput edges and data for graph
    std::cout<<gr<<"\n";
    // Output adjacency matrix for graph
    gr.print_adjacency(std::cout);
  }

  const DistType large=1000000000; // Larger than any possible path length
  DistType min_dist;
  std::vector<bool> zvisited(nvert); // On return will flag whether each vertex has been visited
  std::vector<DistType> distance(nvert); // On return will hold shortest distance for visited vertex
  std::vector<int> parent(nvert); // On return will hold the parent vertex for each visited vertex

  const int nvisited=min_span_tree(gr,
      [](const EdgeDataType &data){return static_cast<DistType>(data);}, large,
      min_dist,zvisited,distance,parent);


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

  if (nvisited==nvert) {
    std::cout<<"Minimum spanning tree path length= "<<min_dist<<"\n";

    int count=0;
    DistType len=0;
    for (int i=0; i<nvert; i++) {
      if (parent[i]>=0) {
        std::cout<<" Edge "<<parent[i]<<"-"<<i<<" Weight "<<distance[i]<<"\n";
        count++;
        len+=distance[i];
      }
    }
    std::cout<<"count="<<count<<" length="<<len<<"\n";
  } else {
    std::cout<<"Graph not full connected\n";
  }

  return 0;
}
