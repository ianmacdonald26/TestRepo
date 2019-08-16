/*
 * graph.cpp
 *
 *  Created on: 13 Aug 2019
 *      Author: Ian
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
// A graph is described by the neighbours onnected to each vertex.
// The class "Neighbour" is used to describe each a single neighbour and the
// data associated with the corresponding edge.
// The edge data type is templated (T)
//***************************************************************************
template<class T> class Neighbour;
template<class T> std::ostream &operator<<(std::ostream &out,
		Neighbour<T> const neigh);

template<class T>
class Neighbour {
public:
	// Constructor
	Neighbour() : vertex(-1), data() {}
	Neighbour(int vertex, T const &data): vertex(vertex), data(data) {};
	// Get vertex number
	int get_vertex() const {return vertex;};
	// Get reference to edge data
	T &get_data() {return data;};
	// Get const reference to edge data
	const T &get_data() const {return data;};

	friend std::ostream &operator<< <>(std::ostream &out,
		Neighbour<T> const neigh);
private:
	int vertex; // Connecting vertex
	T data;     // Edge data
};

template<class T> inline std::ostream &operator<<(std::ostream &out,
		Neighbour<T> const neigh) {
	out<<"("<<neigh.vertex<<" "<<neigh.data<<") ";
	return out;
}

//***************************************************************************
// The class "NeighbourSet" holds the set of neighbours for a particular vertex.
// This implementation stores a std::vector of the "Neighbour" objects.
// Iterators are provided to allow the the user to iterate over the set of neighbours
// The class is templated over edge data type (T)
//***************************************************************************
template<class T> class NeighbourSet;
template<class T> class NeighbourSet_iterator;
template<class T> class NeighbourSet_const_iterator;
template<class T> std::ostream &operator<<(std::ostream &out,
		NeighbourSet<T> const list);

template<class T>
class NeighbourSet {
public:
	using iterator=NeighbourSet_iterator<T>;
	using const_iterator=NeighbourSet_const_iterator<T>;
	// Constructor
  NeighbourSet(int nvert){}; //nvert not needed for this implementation
	// Add neighbour
	void add_neighbour(int vertex, const T &data) {
		neighbours.push_back(Neighbour<T>(vertex,data));
	}
	// Find neighbour and return pointer to data (nullptr if not found)
	const T *find_neighbour(int i) {
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
		NeighbourSet<T> const list);

private:
	// Vector of neighbours
	std::vector<Neighbour<T>> neighbours;
};

template<class T>
inline std::ostream &operator<<(std::ostream &out,
		NeighbourSet<T> const list) {
	for (auto neigh: list.neighbours)
		out<<neigh<<" ";
	return out;
}

// iterator class
template<class T>
class NeighbourSet_iterator {
public:
	using self_type=NeighbourSet_iterator;
	using iterator_category=std::forward_iterator_tag;
	using value_type=Neighbour<T>;
	using difference_type=std::ptrdiff_t;
	using pointer=Neighbour<T>*;
	using reference=Neighbour<T>&;
	// Constructor
	NeighbourSet_iterator(
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
	typename std::vector<Neighbour<T>>::iterator it;
};

// const_iterator class
template<class T>
class NeighbourSet_const_iterator {
public:
	using self_type=NeighbourSet_const_iterator;
	using iterator_category=std::forward_iterator_tag;
	using value_type=Neighbour<T>;
	using difference_type=std::ptrdiff_t;
	using pointer=const Neighbour<T>*;
	using reference=const Neighbour<T>&;
	// Constructor
	NeighbourSet_const_iterator(
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
	typename std::vector<Neighbour<T>>::const_iterator it;
};

//***************************************************************************
// The class "NeighbourSetAdj" holds the set of neighbours for a particular vertex.
// This implementation stores a std::vector with an entry for the all possible neighbours.
// Only the entries corresponding to actual neighbours have valid vertex numbers set.
// Iterators are provides to allow the the user to iterate over the set of neighbours
// The iterators are more complicated and inefficient for this implementation,
// because they need to skip over the entries that are not neighbours.
// The class is templated over edge data type (T).
//***************************************************************************
template<class T> class NeighbourSetAdj;
template<class T> class NeighbourSetAdj_iterator;
template<class T> class NeighbourSetAdj_const_iterator;
template<class T> std::ostream &operator<<(std::ostream &out,
		NeighbourSetAdj<T> const list);

template<class T>
class NeighbourSetAdj {
public:
	using iterator=NeighbourSetAdj_iterator<T>;
	using const_iterator=NeighbourSetAdj_const_iterator<T>;
	// Constructor
	NeighbourSetAdj(int nvert) : neighbours(nvert) {};
	void add_neighbour(int vertex, const T &data) {
		neighbours[vertex]=Neighbour<T>(vertex,data);
	}
	// Find neighbour and return pointer to data
	const T *find_neighbour(int i) {
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
		NeighbourSetAdj<T> const list);

private:
	// Vector of neighbours
	std::vector<Neighbour<T>> neighbours;
};

template<class T>
inline std::ostream &operator<<(std::ostream &out,
		NeighbourSetAdj<T> const set) {
	for (auto neigh: set.neighbours)
		if (neigh.get_vertex()>0) out<<neigh<<" ";
	return out;
}

// iterator class
template<class T>
class NeighbourSetAdj_iterator {
public:
	using self_type=NeighbourSetAdj_iterator;
	using iterator_category=std::forward_iterator_tag;
	using value_type=Neighbour<T>;
	using difference_type=std::ptrdiff_t;
	using pointer=Neighbour<T>*;
	using reference=Neighbour<T>&;
	// Constructor
	NeighbourSetAdj_iterator(
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
	typename std::vector<Neighbour<T>>::iterator it;
	typename std::vector<Neighbour<T>>::iterator itend;
};

// const_iterator class
template<class T>
class NeighbourSetAdj_const_iterator {
public:
	using self_type=NeighbourSetAdj_const_iterator;
	using iterator_category=std::forward_iterator_tag;
	using value_type=Neighbour<T>;
	using difference_type=std::ptrdiff_t;
	using pointer=const Neighbour<T>*;
	using reference=const Neighbour<T>&;
	// Constructor
	NeighbourSetAdj_const_iterator(
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
	typename std::vector<Neighbour<T>>::const_iterator it;
	typename std::vector<Neighbour<T>>::const_iterator itend;
};

//***************************************************************************
//***************************************************************************
template<class T> class Graph;
template<class T> std::ostream &operator<<(std::ostream &out, const Graph<T> &g);

template<class T>
class Graph {
public:
	//Constructor
  Graph(int nvert): nvert(nvert), g(nvert, T(nvert)) {};
  T &operator[](int i){return g[i];};
  const T &operator[](int i) const {return g[i];};
  int get_nvert() const {return nvert;};
  friend std::ostream &operator<< <>(std::ostream &out, const Graph<T> &g);
protected:
	int nvert; // Number of vertices
	std::vector<T> g;
};

template<class T>
inline std::ostream &operator<<(std::ostream &out, const Graph<T> &g){
	for (int i=0; i<g.nvert; i++)
		out<<i<<": "<<g[i]<<"\n";
	return out;
}

//***************************************************************************
//***************************************************************************
template<class T, class E, class D>
class Dijkstra : public Graph<T> {
public:
    Dijkstra(int nvert): Graph<T>(nvert) {};
    int shortest_distance(
        int start_vert,  // Vertex to start from
        int end_vert, // Stop when shortest distance has been found for
                      // this vertex
        // Function to map edge data to edge weight
        std::function<D(const E &)> const &weight,
        D large, // Distance larger than any possible path length
        std::vector<bool> &zvisited, // Set to true when vertex has been visited
        std::vector<D> &distance, // Shortest path distance after node has been visited
        std::vector<int> &parent // Previous vertex for shortest distance path
    ) {
      int n_visited=0;
      for (int i=0; i<this->nvert; i++) {
        zvisited[i]=false;
        distance[i]=large;
        parent[i]=-1;
      }
      D cdist;
      int cvert;
      using pair=std::tuple<D,int>;
      std::priority_queue<pair,std::vector<pair>,std::greater<pair>> pq;
      pq.push(std::make_pair(0.0,start_vert));
      distance[start_vert]=0.0;
      while(!pq.empty()) {
        std::tie(cdist,cvert)=pq.top();
        pq.pop();
        if (!zvisited[cvert]) {
          zvisited[cvert]=true;
          n_visited++;
          if ((n_visited==this->nvert) || (cvert==end_vert)) break;
          for (auto e: this->g[cvert]) {
            const int neigh=e.get_vertex();
            if (!zvisited[neigh]) {
              const D neigh_dist=cdist+weight(e.get_data());
              //std::cout<<cvert<<" "<<neigh<<" "<<distance[neigh]<<
              //" "<<neigh_dist<<"\n";
              if (neigh_dist<distance[neigh]) {
                pq.push(std::make_pair(neigh_dist,neigh));
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
//***************************************************************************
template<class G, class E>
class RandomUndirectedGraph : public G {
public:
	RandomUndirectedGraph(int nvert,
			std::default_random_engine &gen,
			double p,
			std::function<E(std::default_random_engine &)> const &ran_data) 
			: G(nvert) {

			auto dist=std::uniform_real_distribution<double>(0.0,1.0);
		for (int i=0; i<nvert; i++) {
			for (int j=i+1; j<nvert; j++) {
				if (dist(gen)<p) {
					const E data=ran_data(gen);
					this->g[i].add_neighbour(j,data);
					this->g[j].add_neighbour(i,data);
				}
			}
		}
	}
};

} // namespace graph


template <typename T,
          typename Distribution = typename std::conditional<
              std::is_integral<T>::value,
              std::uniform_int_distribution<T>,
              std::uniform_real_distribution<T>>::type>
Distribution get_right_distribution(const T a, const T b) {
    return Distribution(a, b);
}

template<class T> void test(int);

template<class container, class edge_data, class dist_type>
void montecarlo(int nvert, double prob, dist_type wmin, dist_type wmax,
                int nsamples, int outfreq, unsigned seed);

int main() {

	//test<int>(10); // Test routines
	//test<double>(10); // Test routines

	// Use double distances
	using dist_type=double;
	using edge_data=dist_type;

	//Also works with int distances
	//using dist_type=int;
	//using edge_data=dist_type;
	
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
	
	std::cout<<"Edge list storage\n";
	montecarlo<graph::NeighbourSet<edge_data>, edge_data, dist_type>
                  (nvert, prob, wmin, wmax, nsamples, outfreq, seed);
    
  std::cout<<"Adjacency vector storage\n";
  montecarlo<graph::NeighbourSetAdj<edge_data>, edge_data, dist_type>
                  (nvert, prob, wmin, wmax, nsamples, outfreq, seed);

	return 0;
}

template<class container, class edge_data, class dist_type>
void montecarlo(int nvert, double prob, dist_type wmin, dist_type wmax,
                int nsamples, int outfreq, unsigned seed) {

	// Create random generator
	std::default_random_engine gen(seed);

	const dist_type large=1000000000;
	const int start=0; // Start vertex
	const int end=-1;  // End vertex set to -1 so Dijkstra doesn't stop until visited all possible vertices
	double accumavg=0.0; // Running average of shortest distance
	int ntot=0; // Total number of successful (connected) samples

	// Vector data set by shortest_distance method
	std::vector<bool> zvisited(nvert); // On return will flag whether each vertex has been visited
	std::vector<dist_type> distance(nvert); // On return will hold shortest distance for visited vertex
	std::vector<int> parent(nvert); // On return will hold the parent vertex for each visited vertex

	using Dijkstra_graph=graph::Dijkstra<container, edge_data, dist_type>;

	auto tstart = std::chrono::high_resolution_clock::now();
	for (int count=0; count<nsamples; count++) {

		// Create the random Graph
		graph::RandomUndirectedGraph<Dijkstra_graph, edge_data> rg(nvert,  gen, prob,
				[wmin,wmax](std::default_random_engine &gen)
				{
					auto dist=get_right_distribution<dist_type>(wmin,wmax);
					return static_cast<edge_data>(dist(gen));
				}
		);

		// Compute shortest distances
		const int nvisited=rg.shortest_distance(start,end,
				[](const edge_data &data){return static_cast<dist_type>(data);}, large,
				zvisited,distance,parent);

		// If fully connected
		if (nvisited==nvert) {
			// Compute average shortest distance
			const double avgdist=std::accumulate(distance.cbegin()+1,
							distance.cend(), static_cast<dist_type>(0))/(nvert-1.0);
			ntot++;
			accumavg=accumavg*((ntot-1.0)/ntot)+avgdist/ntot;
			if ((outfreq>0) && (ntot%outfreq==0))
				std::cout<<"ntot="<<ntot<<" avg. shortest dist="<<accumavg<<"\n";
		}
	}
	auto tfinish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = tfinish - tstart;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";

	// Print out results
	std::cout<<nsamples-ntot<<" graphs out of "<<nsamples<<" were not fully connected\n";
	std::cout<<"average shortest distance for "<<nvert<<" vertices with probability "<<
			prob<<" and weight in ("<<wmin<<","<<wmax<<") is "<<
			accumavg<<" ("<<ntot<<" samples) \n";
}

template<class T>
void test(int nvert) {

	using dist_type=T;
	using edgedata=dist_type;
	const dist_type large=1000000000;

	unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen(seed);

	using Dijkstra_graph=graph::Dijkstra<graph::NeighbourSet<edgedata>, edgedata, dist_type>;

	graph::RandomUndirectedGraph<Dijkstra_graph, edgedata> rg(nvert, gen, 0.4,
			[](std::default_random_engine &gen)
			{
				auto dist=get_right_distribution<dist_type>(1,10);
				return static_cast<edgedata>(dist(gen));
			}
	);

	std::cout<<rg<<"\n";

	for (int i=0; i<nvert; i++) {
		std::cout<<i<<": ";
		for (int j=0; j<nvert; j++) {
			const edgedata *ed=rg[i].find_neighbour(j);
			if (ed) {
				std::cout<<"("<<j<<" "<<*ed<<")";
			} else {
				std::cout<<"("<<j<<" false)";
			}
		}
		std::cout<<"\n";
	}

	std::vector<bool> zvisited(nvert);
	std::vector<dist_type> distance(nvert);
	std::vector<int> parent(nvert);

	const int nvisited=rg.shortest_distance(0,-1,
			[](const edgedata &data){return static_cast<dist_type>(data);}, large,
			zvisited,distance,parent);

	std::cout<<"nvisited="<<nvisited<<"\n";
	for (int i=1; i<nvert; i++) std::cout<<zvisited[i]<<" ";
	std::cout<<"\n";
	for (int i=1; i<nvert; i++) std::cout<<distance[i]<<" ";
	std::cout<<"\n";
	for (int i=1; i<nvert; i++) std::cout<<parent[i]<<" ";
	std::cout<<"\n";

	if (nvisited==nvert) {
		const double avgdist=std::accumulate(distance.cbegin()+1,
						distance.cend(), static_cast<dist_type>(0))/(nvert-1.0);
		std::cout<<"average shortest distance is "<<avgdist<<"\n";
	} else {
		std::cout<<"not fully connected \n";
	}
}



