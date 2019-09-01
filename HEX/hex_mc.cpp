/*
 * hex_mc.cpp
 *
 *  Created on: 1 Sep 2019
 *      Author: Ian
 */
#include "hex.h"
#include <iostream>
#include <iomanip>
#include <algorithm>  //std::shuffle
#include <chrono>

// Play a computer move
// Returns zero for unfinished, positive for win, negative for exit
int hex::computer_move() {
  // return human_move(); // Play as human

  // Get array of available moves
  std::vector<int> mlist(nrem);
  int cnt=0;
  for(int i=0; i<nloc; i++) {
     if (board[i]==colour::BLANK) mlist[cnt++]=i;
  }

  // Setup permutation array
  std::vector<int> perm(nrem);
  for (int i=0; i<nrem; i++) perm[i]=i;

  // Number of moves that will be colour turn
  const int nfirst= nrem%2==0 ? nrem/2 : nrem/2+1;
  // std::cout<<"nfirst="<<nfirst<<"\n";

  // Array to count number of wins and total number of runs
  std::vector<int> totcount(nrem,0);
  std::vector<int> wincount(nrem,0);

  // Record the start time of simulation
  auto tstart = std::chrono::high_resolution_clock::now();

  const int maxcnt=1000000;
  for (int cnt=0; cnt<maxcnt; cnt++) {

    // Shuffle permutation vector
    std::shuffle(perm.begin(),perm.end(),gen);
    /*for (auto v: mlist) std::cout<<std::setw(3)<<v<<" ";
    std::cout<<"\n";
    for (auto v: perm) std::cout<<std::setw(3)<<v<<" ";
    std::cout<<"\n";*/

    // Fill blank spaces based on permutation vector.
    // Entries 0,...,nfirst-1 are colour turn, the remainder are nextturn
    for (int i=0; i<nfirst; i++) board[mlist[perm[i]]]=turn;
    for (int i=nfirst; i<nrem; i++) board[mlist[perm[i]]]=nextturn;

    // Has colour turn won?
    const bool z=check_for_win(turn,false);

    if (z) {
      // Win
      // Increment totcount and wincount for the colour turn moves
      for (int i=0; i<nfirst; i++) {
        int const loc=perm[i];
        totcount[loc]++;
        wincount[loc]++;
      }
    } else {
      // Loss
      // Increment totcount for the colour turn moves
      for (int i=0; i<nfirst; i++) {
        int const loc=perm[i];
        totcount[loc]++;
      }
    }
  }

  // Record the end time of simulation
  auto tend = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = tend - tstart;
  std::cout << maxcnt<<" simultations. Elapsed time " << elapsed.count() << " s\n";

  // Compute statistics
  // Find move with highest win/run ratio
  int mintot=maxcnt, maxtot=0;
  int minwin=maxcnt, maxwin=0;
  double meantot=0 ,meanwin=0.0;
  double minratio=1.0, maxratio=0.0;
  int best=0;
  for (int i=0; i<nrem; i++) {
    mintot=std::min(mintot,totcount[i]);
    minwin=std::min(minwin,wincount[i]);
    maxtot=std::max(maxtot,totcount[i]);
    maxwin=std::max(maxwin,wincount[i]);
    meantot+=totcount[i];
    meanwin+=wincount[i];
    // Compute win to run ratio
    const double ratio= totcount[i]>0 ?
        static_cast<double>(wincount[i])/totcount[i] : 0.0;
    minratio=std::min(minratio,ratio);
    // Is ration best found so far
    if (ratio>maxratio) {
      maxratio=ratio; // Update best ratio and store move
      best=i;
    }
  }
  meantot=meantot/nrem;
  meanwin=meanwin/nrem;

  // Get column and row of best move
  int row,col;
  std::tie(col,row)=get_pos(mlist[best]);

  // Output statistics
  std::cout<<"mintot="<<mintot<<" maxtot="<<maxtot<<" meantot="<<meantot<<"\n";
  std::cout<<"minwin="<<minwin<<" maxwin="<<maxwin<<" meanwin="<<meanwin<<"\n";
  std::cout<<"minratio="<<minratio<<" maxratio="<<maxratio<<
      " best move: col="<<col<<" row="<<row<<"\n";

  // Set board back to original state
  for (int i=0; i<nrem; i++) board[mlist[i]]=colour::BLANK;

  // Make best move.
  std::cout<<"Computer "<<get_turn()<<"> "<<col<<" "<<row<<"\n";
  // The method move returns true if game has been won
  return move(col,row) ? 1 : 0;
}
