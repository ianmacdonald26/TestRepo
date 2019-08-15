/*
 * graph.cpp
 *
 *  Created on: 13 Aug 2019
 *      Author: ian
 */
#include<iostream>
#include<vector>
#include<random>
#include<functional>
#include<chrono>
#include<queue>
#include<tuple>
#include<numeric>

namespace graph {

//***************************************************************************
template<class T> class neighbour;
template<class T> std::ostream &operator<<(std::ostream &out,
		neighbour<T> const neigh);

template<class T>
class neighbour {
public:
	// Constructor
	neighbour(int vertex, T const &data): vertex(vertex), data(data) {};
	// Get vertex number
	int get_vertex() const {return vertex;};
	// Get reference to edge data
	T &get_data() {return data;};
	// Get const reference to edge data
	const T &get_data() const {return data;};

	friend std::ostream &operator<< <>(std::ostream &out,
		neighbour<T> const neigh);
private:
	int vertex; // Connecting vertex
	T data;     // Edge data
};

template<class T> inline std::ostream &operator<<(std::ostream &out,
		neighbour<T> const neigh) {
	out<<"("<<neigh.vertex<<" "<<neigh.data<<") ";
	return out;
}

//***************************************************************************
template<class T> class neighbour_list;
template<class T> class neighbour_list_iterator;
template<class T> class const_neighbour_list_iterator;
template<class T> std::ostream &operator<<(std::ostream &out,
		neighbour_list<T> const list);

template<class T>
class neighbour_list {
public:
	using iterator=neighbour_list_iterator<T>;
	using const_iterator=const_neighbour_list_iterator<T>;
	// Constructor
	neighbour_list(int nvert){}; //nvert not needed for this implementation
	// Add neighbour
	void add_neighbour(int vertex, const T &data) {
		neighbours.push_back(neighbour<T>(vertex,data));
	}
	// Find neighbour and return pointer to data
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
		neighbour_list<T> const list);

private:
	// Vector of neighbours
	std::vector<neighbour<T>> neighbours;
};

template<class T>
inline std::ostream &operator<<(std::ostream &out,
		neighbour_list<T> const list) {
	for (auto neigh: list.neighbours)
		out<<neigh<<" ";
	return out;
}

// iterator class
template<class T>
class neighbour_list_iterator {
public:
	using self_type=neighbour_list_iterator;
	using iterator_category=std::forward_iterator_tag;
	using value_type=neighbour<T>;
	using difference_type=std::ptrdiff_t;
	using pointer=neighbour<T>*;
	using reference=neighbour<T>&;
	// Constructor
	neighbour_list_iterator(
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
	typename std::vector<neighbour<T>>::iterator it;
};

// const_iterator class
template<class T>
class const_neighbour_list_iterator {
public:
	using self_type=const_neighbour_list_iterator;
	using iterator_category=std::forward_iterator_tag;
	using value_type=neighbour<T>;
	using difference_type=std::ptrdiff_t;
	using pointer=const neighbour<T>*;
	using reference=const neighbour<T>&;
	// Constructor
	const_neighbour_list_iterator(
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
	typename std::vector<neighbour<T>>::const_iterator it;
};

//***************************************************************************
template<class T> class graph;
template<class T> std::ostream &operator<<(std::ostream &out, const graph<T> &g);

template<class T>
class graph {
public:
	//Constructor
	graph(int nvert): nvert(nvert), g(nvert, neighbour_list<T>(nvert)) {};
	neighbour_list<T> &operator[](int i){return g[i];};
	const neighbour_list<T> &operator[](int i) const {return g[i];};
	int get_nvert() const {return nvert;};
	friend std::ostream &operator<< <>(std::ostream &out, const graph<T> &g);
protected:
	int nvert; // Number of vertices
	std::vector<neighbour_list<T>> g;
};

template<class T>
inline std::ostream &operator<<(std::ostream &out, const graph<T> &g){
	for (int i=0; i<g.nvert; i++)
		out<<i<<": "<<g[i]<<"\n";
	return out;
}

//***************************************************************************
template<class T, class D>
class Dijkstra : public graph<T> {
public:
	Dijkstra(int nvert): graph<T>(nvert) {};
	int shortest_distance(
			int start_vert,
			int end_vert,
		    std::function<D(const T &)> const &weight,
			D large,
			std::vector<bool> &zvisited,
			std::vector<D> &distance,
			std::vector<int> &parent) {
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
template<class G, class E>
class random_undirected_graph : public G {
public:
	random_undirected_graph(int nvert,
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

//***************************************************************************

template <typename T,
          typename Distribution = typename std::conditional<
              std::is_integral<T>::value,
              std::uniform_int_distribution<T>,
              std::uniform_real_distribution<T>>::type>
Distribution get_right_distribution(const T a, const T b) {
    return Distribution(a, b);
}

template<class T> void test(int);

int main() {

	test<int>(10); // Test routines
	test<double>(10); // Test routines

	// Use double distances
	using dist_type=double;
	using edgedata=dist_type;

	//Also works with int distances
	//using dist_type=int;
	//sing edgedata=dist_type;
	
	const dist_type large=1000000000;

	int nvert; // Number of vertices
	double prob; // Probability of each possible edge existing
	dist_type wmin; // Minimum edge weight >0
	dist_type wmax; // Maximum edge weigth >=wmin
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
	// Create random generator
	std::default_random_engine gen(seed);

	const int start=0; // Start vertex
	const int end=-1;  // End vertex set to -1 so Dijkstra doesn't stop until visited all possible vertices
	double accumavg=0.0; // Running average of shortest distance
	int ntot=0; // Total number of successful (connected) samples

	// Vector data set by shortest_distance method
	std::vector<bool> zvisited(nvert); // On return will flag whether each vertex has been visited
	std::vector<dist_type> distance(nvert); // On return will hold shortest distance for visited vertex
	std::vector<int> parent(nvert); // On return will hold the parent vertex for each visited vertex

	using Dijkstra_graph=graph::Dijkstra<edgedata, dist_type>;

	for (int count=0; count<nsamples; count++) {

		// Create the random graph
		graph::random_undirected_graph<Dijkstra_graph, edgedata> rg(nvert,  gen, prob,
				[wmin,wmax](std::default_random_engine &gen)
				{
					auto dist=get_right_distribution<dist_type>(wmin,wmax);
					return static_cast<edgedata>(dist(gen));
				}
		);

		// Compute shortest distances
		const int nvisited=rg.shortest_distance(start,end,
				[](const edgedata &data){return static_cast<dist_type>(data);}, large,
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
	// Print out results
	std::cout<<nsamples-ntot<<" graphs out of "<<nsamples<<" were not fully connected\n";
	std::cout<<"average shortest distance for "<<nvert<<" vertices with probability "<<
			prob<<" and weight in ("<<wmin<<","<<wmax<<") is "<<
			accumavg<<" ("<<ntot<<" samples) \n";

	return 0;
}

template<class T>
void test(int nvert) {

	using dist_type=T;
	using edgedata=dist_type;
	const dist_type large=1000000000;

	unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine gen(seed);

	using Dijkstra_graph=graph::Dijkstra<edgedata, dist_type>;

	graph::random_undirected_graph<Dijkstra_graph, edgedata> rg(nvert, gen, 0.4,
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



