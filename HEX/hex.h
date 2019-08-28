#ifndef HEX_H
#define HEX_H

#include <iostream>
#include <vector>
#include <string>

enum colour: char {
  BLANK='.',
  BLUE='b',
  RED='r',
  BLUEWIN='B',
  REDWIN='R'
};

class hex {
  public:
    hex(int dim, bool zfirst=true);
    hex(const std::string &filename);
    int next_move();
    void print_moves(std::ostream &out) const;
  private:
    int dim;
    colour human=colour::BLUE;
    int nloc;
    colour turn=colour::BLUE;
    std::vector<colour> board;
    std::vector<std::pair<int,int>> moves;

    int get_loc(int col, int row) const {return row*dim+col;};
    std::string get_turn() const {return (turn==colour::BLUE) ? "BLUE" : "RED";};
    void move(int col, int row);
    void next_turn() {turn = (turn==colour::BLUE) ? colour::RED : colour::BLUE;};
    void display(std::ostream &out, const std::vector<colour> &board) const;
    int human_move();
    friend std::ostream &operator<<(std::ostream &out, const hex &h);
};

#endif
