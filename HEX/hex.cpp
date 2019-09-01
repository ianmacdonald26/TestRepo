/*
 * hex.cpp
 *
 *  Created on: 28 Aug 2019
 *      Author: Ian
 */
#include "hex.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include "Dijkstra.h"

// Constructor. Sets up board an plays game
hex::hex(
    int dim,      // Dimesnion of board
    bool zfirst,  // true if human is going first
    const std::vector<std::pair<int,int>> &mov  // Vector of pre-existing moves
) :
  dim(dim),
  human(zfirst ? colour::BLUE : colour::RED),
  nloc(dim*dim),
  nrem(nloc),
  board(nloc, colour::BLANK),
  graph(nloc),
  zvisited(nloc),
  dist(nloc),
  parent(nloc),
  zend(nloc)
{
  // Randomize seed for generator
  //unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed=123456;
  gen.seed(seed);

  create_graph(); // Create graph corresponding to board

  // Play pre-existing moves
  for (auto m: mov) {
    const int col=m.first;
    const int row=m.second;
    if ((col<0)||(col>dim-1)||(row<0)||(row>dim)) {
      std::cout<<"Illegal position col="<<col<<" row="<<row<<"\n";
      return;
    }
    if (move(col,row))
    {
      std::cout<<*this<<"\n"; // Display board
      std::cout<<get_turn()<<" has won\n";
      return;
    }
  }

  std::cout<<*this<<"\n"; // Display the board

  play(); // Play game

}

// Method to play a complete game
int hex::play() {
  int status;
  while((status=next_move())==0); // Pay moves until finished

  // Is there a winner or has game been exited
  if (status>0) {
    std::cout<<"******* "<<get_turn()<<" has won ******\n";
  } else {
    std::cout<<"Game Exited\n";
  }

  // Print out moves made
  print_moves(std::cout);

  // Return positive for winner, negative for exited
  return status;
}

// Play a human or computer move depending on turn
int hex::next_move() {
  int status;

  if (turn==human) {
    status=human_move();
  } else {
    status=computer_move();
  }
  // status will be zero for unfinished, positive for win, negative for exit

  if (status>=0) std::cout<<*this<<"\n";

  return status;
}

// Play a human move by asking for row and column
// Returns zero for unfinished, positive for win, negative for exit

int hex::human_move() {
  int row, col;
  int status;
  int loc;
  do { // Read in "col row" until valid move has been made
    status=0;
    std::string line;
    std::cout<<get_turn()<<"> "; // Input prompt
    // Read in next line, ignoring empty lines
    do {
      std::getline(std::cin, line);
    } while (line.empty());
    if ((line=="q")||(line=="Q")) return -1; // Quit game

    if (!(std::stringstream(line)>>col>>row)) { // Try to read col row from string
      std::cout<<"Input error\n";
      status=-1;
    } else {
      if ((row<0)||(row>=dim)) {
        std::cout<<"Invalid row no.\n";
        status=-1;
      }
      if ((col<0)||(col>=dim)) {
        std::cout<<"Invalid column no.\n";
        status=-1;
      }
    }
    // col and row is valid, check location is still empty
    if ((status==0)&&(board[loc=get_loc(col,row)]!=colour::BLANK)) {
        std::cout<<"Location is not empty\n";
        status=-1;
    }
  } while (status<0);

  // Make move.
  // The method move returns true if game has been won
  return move(col,row) ? 1 : 0;
}

// Make (col,row) move on board that has already been check to be valid
bool hex::move(int col, int row) {
  const int loc=get_loc(col,row);
  board[loc]=turn;                          // Update the board
  moves.push_back(std::make_pair(col,row)); // Add to list of moves made
  nrem--; // Decrement number of remaining empty locations
  if (check_for_win(turn,true)) return true; // Check if move has won game
  next_turn();
  return false;
}

// Print an ascii representation of the Hex board
void hex::display(std::ostream &out, const std::vector<colour> &board) const {
  int indent=0;
  for (int row=dim-1; row>0; row--) {
    for (int k=0; k<indent; k++) out<<" ";
    out<<" B"<<std::setw(2)<<row<<" ";
    for (int col=0; col<dim-1; col++) {
      out<<static_cast<char>(board[get_loc(col,row)])<<"  -  ";
    }
    out<<static_cast<char>(board[get_loc(dim-1,row)])<<"\n";
    indent++;
    for (int k=0; k<indent+5; k++) out<<" ";
    for (int col=0; col<dim-1; col++) out<<"\\  /  ";
    out<<"\\\n";
    indent++;
  }
  for (int k=0; k<indent; k++) out<<" ";
  out<<" B"<<std::setw(2)<<0<<" ";
  for (int col=0; col<dim-1; col++) {
    out<<static_cast<char>(board[get_loc(col,0)])<<"  -  ";
  }
  out<<static_cast<char>(board[get_loc(dim-1,0)])<<"\n";
  for (int k=0; k<indent+4; k++) out<<" ";
  for (int col=0; col<dim; col++) out<<std::setw(2)<<col<<"    ";
  out<<"\n";
  for (int k=0; k<indent+5; k++) out<<" ";
  for (int col=0; col<dim; col++) out<<"R     ";
  out<<"\n";
}

// Print the board dimension, whether human is first
// and the list of moves made
void hex::print_moves(std::ostream &out) const {
  out<<dim<<" #dim\n";
  out<<(human==colour::BLUE)<<" #zfirst\n";
  for (int k=0; k<moves.size(); k++) {
    if (k%2==0)
      out<<moves[k].first<<" "<<moves[k].second<<" #Blue\n";
    else
      out<<moves[k].first<<" "<<moves[k].second<<" #Red\n";
  }
}

// Create the graph connectivity
void hex::create_graph() {

  // The six neighbours for interior locations
  for (int col=1; col<dim-1; col++) {
    for (int row=1; row<dim-1; row++) {
      const int n=get_loc(col,row);
      graph[n].add_neighbour(get_loc( col+1, row   ),1);
      graph[n].add_neighbour(get_loc( col-1, row   ),1);
      graph[n].add_neighbour(get_loc( col,   row+1 ),1);
      graph[n].add_neighbour(get_loc( col,   row-1 ),1);
      graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
      graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
    }
  }

  // The four neighbours for bottom row locations (excluding corners)
  for (int col=1; col<dim-1; col++) {
    const int row=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
  }

  // The four neighbours for top row locations (excluding corners)
  for (int col=1; col<dim-1; col++) {
    const int row=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
  }

  // The four neighbours for left edge locations (excluding corners)
  for (int row=1; row<dim-1; row++) {
    const int col=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
  }

  // The four neighbours for right edge locations (excluding corners)
  for (int row=1; row<dim-1; row++) {
    const int col=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
  }

  // The three neighbours for the bottom left corner
  {
    const int col=0;
    const int row=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
  }

  // The two neighbours for the bottom right corner
  {
    const int col=dim-1;
    const int row=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
  }

  // The two neighbours for the top left corner
  {
    const int col=0;
    const int row=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
  }

  // The three neighbours for the top right corner
  {
    const int col=dim-1;
    const int row=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
  }
}

// Checks whether a specific colour has won
// If zpath is true, then the winning path locations are
// changed to BLUWWIN or REDWIN
// Uses Dijkstra shortest path algorithm
bool hex::check_for_win(colour c, bool zpath) {
  colour wc;
  using pair=std::pair<int,int>;
  // Priority queue for Dijkstra
  std::priority_queue<pair, std::vector<pair>, std::greater<pair>> pq;

  const int large=100000;
  // Initialise data used by Dijkstra
  for (int i=0; i<nloc; i++) zend[i]=false;

  if (c==colour::BLUE) {
    // Checking for blue win
    wc=colour::BLUEWIN;
    // Set Dijkstra distance zero vertices to be the blue left edge vertices
    // Set end vertices to be the right edge vertices
    for (int row=0; row<dim; row++) {
      const int n=get_loc(0,row);
      if (board[n]==c) {
        pq.push(std::make_pair(0,n)); // Add vertex to priority q with distance 0
      }
      zend[get_loc(dim-1,row)]=true;
    }
  } else {
    // Checking for red win
    wc=colour::REDWIN;
    // Set Dijkstra distance zero vertices to be the red bottom edge vertices
    // Set end vertices to be the top edge vertices
    for (int col=0; col<dim; col++) {
      const int n=get_loc(col,0);
      if (board[n]==c) {
        pq.push(std::make_pair(0,n)); // Add vertex to priority q with distance 0
      }
      zend[get_loc(col,dim-1)]=true;
    }
  }
  // Apply Dijkstra.
  // Lambda is set up so that egdes are ignored where destination vertex is
  // not colour c, by setting edge length to -1.
  // Routine return the closest end vertex (i.e. one with zend[vert]==true
  // if one is found, otherwise it return negative.
  const int iend=Dijkstra_shortest_distance(graph, pq, zend,
      [this, c](int vert, int &data){return (this->board[vert]==c) ? data : -1;},
      large,
      zvisited, dist, parent);

 /*// Debugging prints
  std::cout<<"loc\n";
  for (int row=dim-1; row>=0; row--) {
    for (int col=0; col<dim; col++) std::cout<<std::setw(3)<<get_loc(col,row)<<" ";
    std::cout<<"\n";
  }
  std::cout<<"zvisited\n";
  for (int row=dim-1; row>=0; row--) {
    for (int col=0; col<dim; col++) std::cout<<std::setw(3)<<zvisited[get_loc(col,row)]<<" ";
    std::cout<<"\n";
  }
  std::cout<<"parent\n";
  for (int row=dim-1; row>=0; row--) {
    for (int col=0; col<dim; col++) std::cout<<std::setw(3)<<parent[get_loc(col,row)]<<" ";
    std::cout<<"\n";
  }
  std::cout<<"dist\n";
  for (int row=dim-1; row>=0; row--) {
    for (int col=0; col<dim; col++) std::cout<<std::setw(3)<<
        dist[get_loc(col,row)]<<" ";
    std::cout<<"\n";
  }*/

  if (zpath&&(iend>=0)) {
    // If Won and zpath, trace back shortest path and change colour to wc
    int  curr=iend;
    while (curr>=0) {
      board[curr]=wc;
      curr=parent[curr]; // Previous vertex along shortest path. -1 indicates start
    }
  }

  // Return true is colour c has won, otherwise false
  return iend>=0;
}

// Print out ascii representation of game using << overload
std::ostream &operator<<(std::ostream &out, const hex &h) {
  h.display(out, h.board);
  return out;
}
