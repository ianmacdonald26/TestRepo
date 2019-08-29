/*
 * Dijkstra.h
 *
 *  Created on: 29 Aug 2019
 *      Author: Ian
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "graph.h"

template<class NeighSetType, class DistType, class EdgeLengthFunc>
int Dijkstra_shortest_distance(
    const graph::Graph<NeighSetType> &graph,
    std::priority_queue<
                        std::pair<int,int>,
                        std::vector<std::pair<int,int>>,
                        std::greater<std::pair<int,int>>
                       > &pq,
    EdgeLengthFunc edge_data_to_distance, // Function to map edge data to edge distance
    // i.e. DistType(int vert, const EdgeDataType &)
    // The following three vectors must be at least length nvertex:
    std::vector<bool> &zvisited,     // Set to true when a vertex has been visited.
    std::vector<DistType> &distance, // Shortest path distance to each vertex (if visited).
    std::vector<int> &parent,        // Previous vertex for shortest distance path (if visited).
    std::vector<bool> &zend
) {
    // Initialise return data
    const int nvert=graph.get_nvert();
    int n_visited=0;

    int cvert; // Current vertex visiting
    DistType cdist; // Shortest distance to current vertex
    while(!pq.empty()) { // Loop until priority queue is empty
      std::tie(cdist,cvert)=pq.top(); // Get lowest distance (top) from priority queue
      pq.pop(); // Remove the top pair from queue
      if (!zvisited[cvert]) { // Ignore if already visited
        zvisited[cvert]=true; // Record as visited
        n_visited++;
        // Finished if all vertices have been visited or current
        // vertex is and end_vertex
        if (zend[cvert]) return cvert;
        if (n_visited==nvert) break;
        // Iterate over the neighbours of current vertex
        for (auto e: graph[cvert]) {
          const int neigh=e.get_vertex();

          if (!zvisited[neigh]) { // Ignore if already visited neighbour
            // Distance to is neighbour is current distance + edge distance
            const DistType edge_len=edge_data_to_distance(neigh,e.get_data());
            // std::cout<<cvert<<" "<<neigh<<" "<<edge_len<<"\n";
            if (edge_len>0) {
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
    return -1;
}



#endif /* DIJKSTRA_H_ */
