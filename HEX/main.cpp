
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

  while (game.next_move()==0);

  std::cout<<game<<"\n";

  std::cout<<"END OF GAME\n";

  game.print_moves(std::cout);

	return 0;
}
