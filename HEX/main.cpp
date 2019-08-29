/*
 * main.cpp
 *
 *  Created on: 28 Aug 2019
 *      Author: Ian
 */
#include "hex.h"
#include <iostream>

int main(int argc, char **argv) {

/*
  int dim;
  do {
    std::cout<<"Enter board size (at least 3) >";
    std::cin>>dim;
  } while (dim<3);

  char ch;
  do {
    std::cout<<"Would you like to go first Y/N? >";
    std::cin>>ch;
  } while ((ch!='y')&&(ch!='Y')&&(ch!='n')&&(ch!='N'));
  const bool zfirst=(ch=='y')||(ch=='Y');

  hex game(dim,zfirst);
*/

  hex game("/home/ian/GIT/TestRepo/HEX/hex.dat");

  std::cout<<game<<"\n";

  int status;
  while (true) {
    status=game.next_move();
    if (status!=0) break;
  }

  if (status>0) std::cout<<game.get_turn()<<" has won\n";

  //game.print_moves(std::cout);

	return 0;
}
