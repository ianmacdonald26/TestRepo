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

// Used to record status of a board location, Value is character for ascii board output
// BLUE/RED also used to record current turn, i.e. turn=colour::RED
enum colour: char {
  BLANK='.',   // Empty position
  BLUE='b',    // Filled by blue
  RED='r',     // Filled by red
  BLUEWIN='B', // Blue winning path
  REDWIN='R'   // Red winning path
};

//***********************************************************************************
// Class to play the game Hex on a dim x dim board.
// constructor sets up board and plays game until there is a winner or the player
// quits (q/Q).
// An optional list of pre-existing moves can be supplied to allow the game
// to restart from a previous unfinished game.
//***********************************************************************************
class hex {
  public:
    hex()=delete; // No default constructor
    // Constructor sets up board and plays game
    hex(
        int dim,               // Dimension hex board
        bool zfirst=true,      // True if human is going first
        // List of previous moves, allowing game to re-start from a previous game
        const std::vector<std::pair<int,int>> &mov=std::vector<std::pair<int,int>>()
        );
  private:
    const int dim;                           // Dimension of hex board
    const colour human;                      // Colour of human player
    const int nloc;                          // Number of locations in board (ndim*dmim)
    int nrem;                                // Number of empty locations
    colour turn=colour::BLUE;                // Colour of current turn, blue is always first
    std::vector<colour> board;               // Current status of each location
    std::vector<std::pair<int,int>> moves;   // Stores previous moves
    graph::Graph<graph::NeighbourSetAdjList<int>> graph;  // Graph of hex board
    // The following four vectors are used by the Dijkstra shortest path routine
    std::vector<bool> zvisited;
    std::vector<int> dist;
    std::vector<int> parent;
    std::vector<bool> zend;

    int get_loc(int col, int row) const {return row*dim+col;}; // Convert column, row to location
    void next_turn() {turn = (turn==colour::BLUE) ? colour::RED : colour::BLUE;}; // Change turn
    std::string get_turn() const {return (turn==colour::BLUE) ? "BLUE" : "RED";}; // Convert turn to string

    int play();                      // Play game
    int next_move();                 // Play next move
    int human_move();                // Play human move
    int computer_move();             // Play computer move
    bool move(int col, int row);     // Play (col, row) move and check for win
    void display(std::ostream &out, const std::vector<colour> &board) const;  // Display ascii board representation
    void print_moves(std::ostream &out) const;  // Print out board size and list of moves
    void create_graph();             // Create graph of hex board
    bool check_for_win(colour c, bool zpath=false); // Check if colour c has won

    friend std::ostream &operator<<(std::ostream &out, const hex &h);
};

#endif /* HEX_H_ */
