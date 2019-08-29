/*
 * hex.cpp
 *
 *  Created on: 28 Aug 2019
 *      Author: Ian
 */
#include "hex.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "Dijkstra.h"

hex::hex(int dim, bool zfirst) :
  dim(dim),
  nloc(dim*dim),
  board(nloc, colour::BLANK),
  human(zfirst ? colour::BLUE : colour::RED),
  graph(nloc),
  zvisited(nloc),
  dist(nloc),
  parent(nloc),
  zend(nloc)
{
  create_graph();
}

hex::hex(std::string filename) {
  std::string line;

  std::cout<<"Reading input from file: "<<filename<<"\n";
  std::ifstream in(filename);
  if (!in) {
    std::cout<<"Unable to open file: "<<filename<<"\n";
    exit(EXIT_FAILURE);
  }

  {
    std::getline(in, line);
    std::stringstream ss(line);
    ss>>dim;
    if (!ss) {
      std::cout<<"Unable to read dimension from file: "<<filename<<"\n";
      exit(EXIT_FAILURE);
    }
  }
  std::cout<<"dim="<<dim<<"\n";
  nloc=dim*dim;
  board.resize(nloc, colour::BLANK);

  graph.resize(nloc);
  create_graph();
  zvisited.resize(nloc);
  dist.resize(nloc);
  parent.resize(nloc);
  zend.resize(nloc);

  {
    bool zfirst;
    std::getline(in, line);
    std::stringstream ss(line);
    ss>>zfirst;
    if (!ss) {
      std::cout<<"Unable to read zfirst from file: "<<filename<<"\n";
      exit(EXIT_FAILURE);
    }
    std::cout<<"zfirst="<<zfirst<<"\n";
    human = zfirst ? colour::BLUE : colour::RED;
  }

  while (true) {
    int row,col;
    std::getline(in, line);
    std::stringstream ss(line);
    ss>>col>>row;
    if (!ss) break;
    std::cout<<"col="<<col<<" row="<<row<<" #"<<get_turn()<<"\n";
    move(col,row);
  }

  in.close();
}

bool hex::check_for_win(colour c, bool zpath) {
  colour *pb=&board[0];
  colour wc;
  std::priority_queue<
                      std::pair<int,int>,
                      std::vector<std::pair<int,int>>,
                      std::greater<std::pair<int,int>>
                     > pq;

  const int large=100000;
  for (int i=0; i<nloc; i++) {
    zvisited[i]=false;
    dist[i]=large;
    parent[i]=-1;
    zend[i]=false;
  }

  if (c==colour::BLUE) {
    wc=colour::BLUEWIN;
    for (int row=0; row<dim; row++) {
      const int n=get_loc(0,row);
      if (board[n]==c) {
        pq.push(std::make_pair(0,n));
        dist[n]=0;
      }
      zend[get_loc(dim-1,row)]=true;
    }
  } else {
    wc=colour::REDWIN;
    for (int col=0; col<dim; col++) {
      const int n=get_loc(col,0);
      if (board[n]==c) {
        pq.push(std::make_pair(0,n));
        dist[n]=0;
      }
      zend[get_loc(col,dim-1)]=true;
    }
  }
  const int iend=Dijkstra_shortest_distance(graph, pq,
      [pb, c](int vert, int &data){return (pb[vert]==c) ? data : -1;},
      zvisited, dist, parent, zend);
/*
  for (int row=dim-1; row>=0; row--) {
    for (int col=0; col<dim; col++) std::cout<<std::setw(3)<<get_loc(col,row)<<" ";
    std::cout<<"\n";
  }
  for (int row=dim-1; row>=0; row--) {
    for (int col=0; col<dim; col++) std::cout<<std::setw(3)<<zvisited[get_loc(col,row)]<<" ";
    std::cout<<"\n";
  }
*/
  if (zpath&&(iend>=0)) {
    int  curr=iend;
    while (curr>=0) {
      board[curr]=wc;
      curr=parent[curr];
    }
  }

  return iend>=0;
}

void hex::create_graph() {
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
  for (int col=1; col<dim-1; col++) {
    const int row=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
  }
  for (int col=1; col<dim-1; col++) {
    const int row=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
  }
  for (int row=1; row<dim-1; row++) {
    const int col=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
  }
  for (int row=1; row<dim-1; row++) {
    const int col=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
  }
  {
    const int col=0;
    const int row=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
    graph[n].add_neighbour(get_loc( col+1, row+1 ),1);
  }
  {
    const int col=dim-1;
    const int row=0;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row+1 ),1);
  }
  {
    const int col=0;
    const int row=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col+1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
  }
  {
    const int col=dim-1;
    const int row=dim-1;
    const int n=get_loc(col,row);
    graph[n].add_neighbour(get_loc( col-1, row   ),1);
    graph[n].add_neighbour(get_loc( col,   row-1 ),1);
    graph[n].add_neighbour(get_loc( col-1, row-1 ),1);
  }
}

bool hex::move(int col, int row) {
  const int loc=get_loc(row,col);
  board[loc]=turn;
  moves.push_back(std::make_pair(col,row));
  if (check_for_win(turn,true)) return true;
  next_turn();
  return false;
}

int hex::next_move() {
  int status;

  if (turn==human) {
    status=human_move();
  } else {
    status=human_move();
  }

  if (status>=0) std::cout<<*this<<"\n";

  return status;
}

int hex::human_move() {
  int row, col;
  int status;
  int loc;
  do {
    status=0;
    std::string line;
    std::cout<<get_turn()<<"> ";
    std::getline(std::cin, line);
    std::stringstream ss(line);
    ss>>col;
    if (ss && (col<0)) return -1;
    ss>>row;
    if (!ss) {
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
    if ((status==0)&&(board[loc=get_loc(col,row)]!=colour::BLANK)) {
        std::cout<<"Location is not empty\n";
        status=-1;
    }
  } while (status<0);

  if (move(row,col)) status=1;

  return status;
}

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

std::ostream &operator<<(std::ostream &out, const hex &h) {
  h.display(out, h.board);
  return out;
}
