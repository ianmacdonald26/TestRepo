/*
 * Prim.cpp
 *
 *  Created on: 21 Aug 2019
 *      Author: Ian MacDonald
 *      
 *      
 */

/* CODE OUTPUT
Reading graph from file: SampleTestData_mst_data.dat
nvert=20
Adjacency List:
0: (1 17)  (2 2)  (3 9)  (4 24)  (5 28)  (6 29)  (7 14)  (8 28)  (9 13)  (10 23)  (11 10)  (12 15)  (13 23)  (14 15)  (15 18)  (16 11)  (17 4)  (18 27)  (19 5)
1: (0 17)  (2 9)  (3 3)  (5 14)  (6 1)  (8 27)  (9 20)  (10 16)  (11 24)  (12 29)  (13 6)  (15 15)  (16 20)  (17 1)  (18 11)  (19 9)
2: (0 2)  (1 9)  (3 21)  (4 21)  (5 29)  (6 13)  (7 19)  (8 16)  (9 1)  (10 11)  (11 4)  (12 12)  (14 26)  (15 5)  (16 25)  (17 12)  (18 5)  (19 24)
3: (0 9)  (1 3)  (2 21)  (4 11)  (6 22)  (7 22)  (8 12)  (9 16)  (11 22)  (12 1)  (13 12)  (15 14)  (16 15)  (17 23)  (18 27)  (19 28)
4: (0 24)  (2 21)  (3 11)  (5 25)  (7 1)  (8 1)  (9 5)  (11 24)  (12 29)  (13 9)  (14 4)  (15 2)  (16 5)  (18 10)  (19 10)
5: (0 28)  (1 14)  (2 29)  (4 25)  (6 1)  (7 17)  (8 22)  (9 7)  (10 20)  (11 7)  (12 22)  (13 16)  (14 11)  (15 22)  (16 2)  (17 23)  (18 1)  (19 20)
6: (0 29)  (1 1)  (2 13)  (3 22)  (5 1)  (8 18)  (9 7)  (11 4)  (12 18)  (13 11)  (14 14)  (15 5)  (16 24)  (17 5)  (18 13)
7: (0 14)  (2 19)  (3 22)  (4 1)  (5 17)  (8 27)  (9 7)  (10 2)  (11 5)  (13 29)  (14 16)  (15 25)  (16 8)  (17 19)  (18 26)  (19 23)
8: (0 28)  (1 27)  (2 16)  (3 12)  (4 1)  (5 22)  (6 18)  (7 27)  (9 3)  (10 3)  (11 26)  (12 9)  (13 25)  (14 16)  (15 7)  (16 4)  (17 23)  (18 7)
9: (0 13)  (1 20)  (2 1)  (3 16)  (4 5)  (5 7)  (6 7)  (7 7)  (8 3)  (11 23)  (12 3)  (13 3)  (14 28)  (15 24)  (16 12)  (17 20)  (18 25)  (19 25)
10: (0 23)  (1 16)  (2 11)  (5 20)  (7 2)  (8 3)  (12 27)  (13 13)  (14 25)  (15 2)  (16 3)  (17 4)  (18 4)  (19 15)
11: (0 10)  (1 24)  (2 4)  (3 22)  (4 24)  (5 7)  (6 4)  (7 5)  (8 26)  (9 23)  (12 1)  (14 1)  (15 20)  (16 20)  (17 22)  (18 19)  (19 28)
12: (0 15)  (1 29)  (2 12)  (3 1)  (4 29)  (5 22)  (6 18)  (8 9)  (9 3)  (10 27)  (11 1)  (13 23)  (14 6)  (15 9)  (16 28)  (17 1)  (18 6)  (19 13)
13: (0 23)  (1 6)  (3 12)  (4 9)  (5 16)  (6 11)  (7 29)  (8 25)  (9 3)  (10 13)  (12 23)  (14 5)  (15 19)  (16 18)  (17 4)  (18 16)  (19 12)
14: (0 15)  (2 26)  (4 4)  (5 11)  (6 14)  (7 16)  (8 16)  (9 28)  (10 25)  (11 1)  (12 6)  (13 5)  (15 6)  (16 27)  (17 15)  (18 1)  (19 28)
15: (0 18)  (1 15)  (2 5)  (3 14)  (4 2)  (5 22)  (6 5)  (7 25)  (8 7)  (9 24)  (10 2)  (11 20)  (12 9)  (13 19)  (14 6)  (16 23)  (17 21)  (18 28)  (19 2)
16: (0 11)  (1 20)  (2 25)  (3 15)  (4 5)  (5 2)  (6 24)  (7 8)  (8 4)  (9 12)  (10 3)  (11 20)  (12 28)  (13 18)  (14 27)  (15 23)  (17 9)  (18 11)  (19 12)
17: (0 4)  (1 1)  (2 12)  (3 23)  (5 23)  (6 5)  (7 19)  (8 23)  (9 20)  (10 4)  (11 22)  (12 1)  (13 4)  (14 15)  (15 21)  (16 9)  (18 20)  (19 9)
18: (0 27)  (1 11)  (2 5)  (3 27)  (4 10)  (5 1)  (6 13)  (7 26)  (8 7)  (9 25)  (10 4)  (11 19)  (12 6)  (13 16)  (14 1)  (15 28)  (16 11)  (17 20)  (19 11)
19: (0 5)  (1 9)  (2 24)  (3 28)  (4 10)  (5 20)  (7 23)  (9 25)  (10 15)  (11 28)  (12 13)  (13 12)  (14 28)  (15 2)  (16 12)  (17 9)  (18 11)

Minimum spanning tree length= 30
 Edge 17-1 Weight 1
 Edge 0-2 Weight 2
 Edge 12-3 Weight 1
 Edge 8-4 Weight 1
 Edge 6-5 Weight 1
 Edge 1-6 Weight 1
 Edge 4-7 Weight 1
 Edge 9-8 Weight 3
 Edge 2-9 Weight 1
 Edge 7-10 Weight 2
 Edge 12-11 Weight 1
 Edge 9-12 Weight 3
 Edge 9-13 Weight 3
 Edge 11-14 Weight 1
 Edge 4-15 Weight 2
 Edge 5-16 Weight 2
 Edge 12-17 Weight 1
 Edge 14-18 Weight 1
 Edge 15-19 Weight 2
count=19 length=30
*/

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
    Graph(int nvert) : nvert(nvert), g(nvert, NeighSetType(nvert)) {};

    // Constructor to read in graph from given file name
    Graph(std::string filename) : nvert(0), g(0, NeighSetType(0)){
      std::cout<<"Reading graph from file: "<<filename<<"\n";
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
// Main routine reads creates a graph using a
// The Prim_min_span_tree routine is called for this graph and the
// Minimum spanning tree (MST) information is written out.
//***************************************************************************
int main() {
  const bool zverbose=false;

  // Use int distances
  using DistType=int;
  using EdgeDataType=DistType; // The edge data consists of the same type

  // Create graph from given file name.
  graph::Graph<graph::NeighbourSetAdjList<EdgeDataType>> gr("SampleTestData_mst_data.dat");
  const int nvert=gr.get_nvert();

  // Write out graph
  std::cout<<"nvert="<<nvert<<"\n";
  std::cout<<"Adjacency List:\n"<<gr<<"\n"; // Ouput edges and data for graph
  if (zverbose) {
    std::cout<<"Adjacency Matrix:\n";
    gr.print_adjacency_matrix(std::cout); // Output adjacency matrix for graph
  }

  const DistType len_large=1000000000; // Larger than any possible path length
  DistType mst_len;
  std::vector<int>      parent(nvert); // (edges (i,parent[i]) i=1,i...,nvert-1
                                       // will be the minimum spanning tree
  std::vector<DistType> length(nvert); // length[i] will be the length of edge (i,parent[i])
  std::vector<bool>     zclosed(nvert);// Stores closed set membership

  // Compute MST using Prim's algorithm
  const int nclosed=Prim_min_span_tree(gr,
      [](const EdgeDataType &data){return static_cast<DistType>(data);}, len_large,
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
      std::cout<<" Edge "<<parent[i]<<"-"<<i<<" Weight "<<length[i]<<"\n";
      count++;
      len+=length[i];
    }
    std::cout<<"count="<<count<<" length="<<len<<"\n";
  } else {
    std::cout<<"Graph not full connected\n";
  }

  return 0;
}
