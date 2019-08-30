/*
 * hex.h
 *
 *  Created on: 28 Aug 2019
 *      Author: Ian
 */
#ifndef HEX_H_
#define HEX_H_

#include <iostream>
#include <vector>
#include <string>
#include "graph.h"

enum colour: char {
  BLANK='.',
  BLUE='b',
  RED='r',
  BLUEWIN='B',
  REDWIN='R'
};

class hex {
  public:
    hex(
        int dim,
        bool zfirst=true,
        const std::vector<std::pair<int,int>> &mov=std::vector<std::pair<int,int>>()
        );
    int next_move();
    void print_moves(std::ostream &out) const;
    std::string get_turn() const {return (turn==colour::BLUE) ? "BLUE" : "RED";};
  private:
    const int dim;
    const colour human=colour::BLUE;
    const int nloc;
    colour turn=colour::BLUE;
    std::vector<colour> board;
    std::vector<std::pair<int,int>> moves;
    std::vector<bool> zvisited;
    std::vector<int> dist;
    std::vector<int> parent;
    std::vector<bool> zend;

    graph::Graph<graph::NeighbourSetAdjList<int>> graph;

    int get_loc(int col, int row) const {return row*dim+col;};
    bool move(int col, int row);
    void next_turn() {turn = (turn==colour::BLUE) ? colour::RED : colour::BLUE;};
    void display(std::ostream &out, const std::vector<colour> &board) const;
    int human_move();
    void create_graph();
    bool check_for_win(colour c, bool zpath=false);
    friend std::ostream &operator<<(std::ostream &out, const hex &h);
};

#endif /* HEX_H_ */
