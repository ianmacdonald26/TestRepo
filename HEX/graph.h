/*
 * graph.h
 *
 *  Created on: 29 Aug 2019
 *      Author: Ian
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>
#include <tuple>
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
    int vertex;        // Connecting vertex
    EdgeDataType data; // Edge data
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
    // Add neighbour, getting data from supplied input stream
    bool add_neighbour(int vertex, std::istream & in) {
      EdgeDataType data;
      const bool z=in>>data; // Read in edge data
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
// The "Graph<NeighSetType>" class is the container of the neighbour sets for
// all vertices.
// The [] operator is overloaded to return a reference to the neighbour set
// for the requested vertex.
//***************************************************************************
template<class NeighSetType> class Graph;
template<class NeighSetType> std::ostream &operator<<(std::ostream &out, const Graph<NeighSetType> &g);

template<class NeighSetType>
class Graph {
  public:
    //Constructor
    Graph(int nvert=0) : nvert(nvert), g(nvert, NeighSetType(nvert)) {};

    // Constructor to read in graph from given file name
    Graph(std::string filename) : nvert(0), g(0, NeighSetType(0)){
      std::cout<<"Reading graph from file: "<<filename<<"\n";
      std::ifstream in(filename);
      in>>nvert;
      resize(nvert);
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

    void resize(int nvert) {
      g.resize(nvert,NeighSetType(nvert));
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

} // namespace graph
#endif /* GRAPH_H_ */
