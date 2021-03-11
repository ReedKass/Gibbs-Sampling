/* Name: Reed Kass-Mullet
 * Date: 2/17/21
 * Class: COMP167
 * Assignment: HW2 - Part 2
 * Description: Gibbs Sampling
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
using namespace std;

struct gibbs_data {
  int s_star;
  vector<int> positions;
};

vector<string> initialize();
void gibbs(vector<string> sequences, int motif_len, int seq_len, int num_seq);
int gibbs_iterate(vector<string> sequences, struct gibbs_data gd, int motif_len);
vector<vector<float>> score_pssm(vector<vector<float>> pssm, vector<string> sequences, struct gibbs_data gd);
int score_s_star(vector<vector<float>> pssm, vector<string> sequences, struct gibbs_data gd);
void print_seqs(struct gibbs_data gd, vector<string> sequences, int motif_len, int seq_len, int num_seq);

int main(int argc, char* argv[]) {
  (void) argc;
  (void) argv;

  vector<string> sequences = initialize();

  int motif_len = 6;
  int seq_len = sequences.at(0).size();
  int num_seq = sequences.size();
  // ======================================
  for(int i = 0; i < 5; i++){
    cout << sequences.at(i) << endl;
  }
  cout << "seq length: " << seq_len << endl;
  cout << "num seq: " << num_seq << endl;
  cout << "motif len:" << motif_len << endl << "===========================\n";
  //=======================================

  gibbs(sequences, motif_len, seq_len, num_seq);

  return 0;
}

/* This function reads in the sequences, it is not designed for software capability,
   and the number of sequences being read in must be chosen by hand via the variable
   controlling the for loop and the size of the vector
   */
vector<string> initialize() {
  vector<string> sequences(5);
  for (int i = 0; i < 5; i++) {
    string temp = "";
    getline(cin, temp);
    getline(cin, temp);
    sequences.at(i) = temp;
  }

  return sequences;
}

/* This function runs the gibbs algorithm. It uses a for loop to follow the iterative
   nature of the algorithm. A full description of this function and how it works is given
   in the README file.
   */
void gibbs(vector<string> sequences, int motif_len, int seq_len, int num_seq) {
  struct gibbs_data gd;

  //We start by selelcting random s* and initial positions
  srand(time(NULL));
  gd.s_star = rand() % num_seq;
  vector<int> positions(num_seq);
  vector<int> prev_positions(num_seq);
  gd.positions = positions;

  //We randomize initial positions and initialize the previous position vector to track
  //the convergence condition
  for (int i = 0; i < num_seq; i++) {
    gd.positions.at(i) = rand() % (seq_len - motif_len);
    prev_positions.at(i) = -1;
  }

  bool matching_s = false;
  while(!matching_s) {
    //Create and analyze PSSM to find best prefix for s*
    int best_prefix = gibbs_iterate(sequences, gd, motif_len);
    cout << "best_prefix " << best_prefix << endl;
    cout << "s* " << gd.s_star << endl; 

    gd.positions.at(gd.s_star) = best_prefix;
    if(prev_positions.at(gd.s_star) == gd.positions.at(gd.s_star)){
      matching_s = true;
    }
    prev_positions.at(gd.s_star) = gd.positions.at(gd.s_star);
    int next_s_star = rand() % num_seq;
    while(gd.s_star == next_s_star) {
      next_s_star = rand() % num_seq;
    }
    gd.s_star = next_s_star;
  }

  print_seqs(gd, sequences, motif_len, seq_len, num_seq);
}

//This functions prints the final motifs chosen from the convergence condition
void print_seqs(struct gibbs_data gd, vector<string> sequences, int motif_len, int seq_len, int num_seq) {
  cout << "ending print\n";
  for(int i = 0; i < num_seq; i++){
    for(int j = 0; j < seq_len; j++){
      if(j == gd.positions.at(i)){
        cout << " [";
      }
      cout << sequences.at(i).at(j);
      if(j == (gd.positions.at(i) + motif_len)) {
        cout << "] ";
      }
    }
    cout << endl;
  }
}

//This function creates the pssm, calls a function to score the pssm, and then
//passes the pssm into a function which finds the best prefix for s*
int gibbs_iterate(vector<string> sequences, struct gibbs_data gd, int motif_len) {

  vector<vector<float>> pssm(4);
  for(int i = 0; i < 4; i++){
    vector<float> freq(motif_len);
    pssm.at(i) = freq;
  }
  pssm = score_pssm(pssm, sequences, gd);

  //We use pssm to find the best fitting motif on s*
  int best_prefix = score_s_star(pssm, sequences, gd);
  return best_prefix;
}

//This function takes the empty pssm structures, creates a frequency matrix in it,
//and then turns that frequency matrix into the pssm which is returned.
vector<vector<float>> score_pssm(vector<vector<float>> pssm, vector<string> sequences, struct gibbs_data gd) {
  int motif_len = pssm.at(0).size();
  int num_seq = sequences.size();
  cout << "ml " << motif_len << " num seq: " << num_seq << endl;
  for(int i = 0; i < motif_len; i ++) {
    for(int j = 0; j < num_seq; j++) {
      if(j != gd.s_star) {
        if(sequences.at(j).at(i + gd.positions.at(j)) == 'A' ) {
          pssm.at(0).at(i) = 1 + pssm.at(0).at(i);
        }
        else if(sequences.at(j).at(i + gd.positions.at(j)) == 'C' ) {
          pssm.at(1).at(i) = 1 + pssm.at(1).at(i);
        }
        else if(sequences.at(j).at(i + gd.positions.at(j)) == 'G' ) {
          pssm.at(2).at(i) = 1 + pssm.at(2).at(i);
        }
        else if(sequences.at(j).at(i + gd.positions.at(j)) == 'T' ) {
          pssm.at(3).at(i) = 1 + pssm.at(3).at(i);
        }
      }
    }
  }

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < motif_len; j++) {
      pssm.at(i).at(j) = 1 + pssm.at(i).at(j);
      pssm.at(i).at(j) = (pssm.at(i).at(j)/(num_seq - 1))/.25;
      cout << pssm.at(i).at(j) << " ";
    }
    cout << endl;
  }
  return pssm;
}

//This function takes the scored pssm and uses it to score and track all potential
//prefix positions in the s* sequence. It returns the integer value representing
//the optimal location for the position of s* to be updated to.
int score_s_star(vector<vector<float>> pssm, vector<string> sequences, struct gibbs_data gd) {
  int motif_len = pssm.at(0).size();
  int max_score;
  int max_prefix;
  for(int i = 0; i < (int) sequences.at(gd.s_star).size() - motif_len; i++) {
    int score;
    for(int j = 0; j < motif_len; j++) {
      int value;
      if(sequences.at(gd.s_star).at(j + i) == 'A') {
        value = pssm.at(0).at(j);
      }
      else if(sequences.at(gd.s_star).at(j + i) == 'C') {
        value = pssm.at(1).at(j);
      }
      else if(sequences.at(gd.s_star).at(j + i) == 'G') {
        value = pssm.at(2).at(j);
      }
      else if(sequences.at(gd.s_star).at(j + i) == 'T') {
        value = pssm.at(3).at(j);
      }
      if (j == 0) {
        score = value;
      } else {
        score = score * value;
      }
    }
    //cout << "score for " << i << " " << score << endl;
    if (i == 0) {
      max_score = score;
      max_prefix = 0;
    } else if(score > max_score) {
      max_score = score;
      max_prefix = i;
    }
  }
  cout << "best score " << max_score << " best prefix " << max_prefix << endl;
  return max_prefix;
}