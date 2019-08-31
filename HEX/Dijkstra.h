/*
 * Dijkstra.h
 *
 *  Created on: 29 Aug 2019
 *      Author: Ian
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "graph.h"

//*************************************************************************************
// Custom version of Dijkstra shortest path algorithm for use with the Hex board game.
// User must supply a std::priority_queue containing the (distance, vertex) start pairs.
// User must supply the array zend, which is set to true for only the terminal vertices.
// Search will terminate if the current vertex is a "terminal vertex" (zend true)
// and will return this vertex number.
// It will return a negative value if no terminal vertex can be reached from the start
// points.
// The distance array will return the shortest distance to each VISITED vertex.
// The parent array will then give the the previous vertex in a shortest distance path,
// thus allowing a full shortest path to be traced back.
//*************************************************************************************
template<class NeighSetType, class DistType, class EdgeLengthFunc>
int Dijkstra_shortest_distance(
    const graph::Graph<NeighSetType> &graph,                   // The graph object;
    std::priority_queue<                                       // Priority queue of starting
                        std::pair<DistType,int>,               // (distance, vertex) pairs
                        std::vector<std::pair<DistType,int>>,
                        std::greater<std::pair<DistType,int>>
                       > &pq,
    const std::vector<bool> &zend,        // Dijkstra will terminate and return the vertex no. if it
                                          // visits a vertex where zend is true.
    EdgeLengthFunc edge_data_to_distance, // Function to map vertex no. and edge data to edge distance
                                          // i.e. DistType(int vert, const EdgeDataType &)
    const DistType large,                 // Value greater than largest possible path length
    // The following three vectors must be at least length nvertex:
    std::vector<bool> &zvisited,     // Set to true when a vertex has been visited.
    std::vector<DistType> &distance, // Shortest path distance to each vertex (if visited).
    std::vector<int> &parent         // Previous vertex for shortest distance path (if visited).
) {
    const int nvert=graph.get_nvert();
    // Initialise returned data
    for (int i=0; i<nvert; i++) {
      zvisited[i]=false;
      distance[i]=large;
      parent[i]=-1;
    }
    int n_visited=0;
    int cvert;           // Current vertex visiting.
    DistType cdist;      // Shortest distance to current vertex.
    while(!pq.empty()) { // Loop until priority queue is empty.
      std::tie(cdist,cvert)=pq.top(); // Get lowest distance (top) from priority queue.
      pq.pop();                       // Remove the top pair from queue.
      if (!zvisited[cvert]) {         // Ignore if already visited.
        zvisited[cvert]=true;         // Record as visited.
        distance[cvert]=cdist;
        n_visited++;
        // Finished if all vertices have been visited or current
        // vertex is and end_vertex
        if (zend[cvert]) return cvert;
        if (n_visited==nvert) break;
        // Iterate over the neighbours of current vertex (only considers edges with len>0)
        for (auto e: graph[cvert]) {
          const int neigh=e.get_vertex();
          if (!zvisited[neigh]) { // Ignore if already visited neighbour
            // Get edge length
            const DistType edge_len=edge_data_to_distance(neigh,e.get_data());
            if (edge_len>0) { // Ignore egde unless the length is positive
              // Distance to neighbour is current distance + edge length
              const DistType neigh_dist=cdist+edge_len;
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
    }
    return -1; // Have not reached a terminal vertex
}
#endif /* DIJKSTRA_H_ */
