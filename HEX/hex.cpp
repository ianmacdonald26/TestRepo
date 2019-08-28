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

hex::hex(int dim, bool zfirst) :
  dim(dim),
  nloc(dim*dim),
  board(nloc, colour::BLANK),
  human(zfirst ? colour::BLUE : colour::RED)
{}

hex::hex(const std::string &filename) {
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
    std::cout<<"dim="<<dim<<"\n";
    nloc=dim*dim;
    board.resize(nloc, colour::BLANK);
  }

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

void hex::move(int col, int row) {
  const int loc=get_loc(row,col);
  board[loc]=turn;
  moves.push_back(std::make_pair(col,row));
  next_turn();
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

  move(row,col);

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
