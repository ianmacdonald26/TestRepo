/*
 * main.cpp
 *
 *  Created on: 28 Aug 2019
 *      Author: Ian
 */
#include "hex.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

int read_game(std::string filename,
              int &dim,
              bool &zfirst,
              std::vector<std::pair<int,int>> &moves);

int main(int argc, char **argv) {

  int dim;
  bool zfirst;
  std::vector<std::pair<int,int>> moves;

  if (argc>1) {

    if (read_game(argv[1], dim, zfirst, moves))
        exit(EXIT_FAILURE);

  } else {

    do {
      std::cout<<"Enter board size (at least 3) >";
      std::cin>>dim;
    } while (dim<3);

    char ch;
    do {
      std::cout<<"Would you like to go first Y/N? >";
      std::cin>>ch;
    } while ((ch!='y')&&(ch!='Y')&&(ch!='n')&&(ch!='N'));
    zfirst=(ch=='y')||(ch=='Y');

  }

  hex game(dim,zfirst,moves);

	return 0;
}

int read_game(std::string filename,
              int &dim,
              bool &zfirst,
              std::vector<std::pair<int,int>> &moves) {

  std::string line;

  std::cout<<"Reading input from file: "<<filename<<"\n";

  std::ifstream in(filename);
  if (!in) {
    std::cout<<"Unable to open file: "<<filename<<"\n";
    return -1;;
  }

  {
    std::getline(in, line);
    std::stringstream ss(line);
    ss>>dim;
    if (!ss) {
      std::cout<<"Unable to read dimension from file: "<<filename<<"\n";
      return -2;
    }
  }
  std::cout<<"dim="<<dim<<"\n";

  {
    std::getline(in, line);
    std::stringstream ss(line);
    ss>>zfirst;
    if (!ss) {
      std::cout<<"Unable to read zfirst from file: "<<filename<<"\n";
      return -3;
    }
    std::cout<<"zfirst="<<zfirst<<"\n";
  }

  while (true) {
    int row,col;
    std::getline(in, line);
    std::stringstream ss(line);
    ss>>col>>row;
    if (!ss) break;
    std::cout<<"col="<<col<<" row="<<row<<"\n";
    moves.push_back(std::make_pair(col,row));
  }

  in.close();
  return 0;
}
