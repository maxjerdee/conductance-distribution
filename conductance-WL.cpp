// This script uses flat histrogram methods to estimate the distribution of conductances among
// possible partitions of a generic graph. 
// The graph may be input as an edgelist in edges.txt
// The (logarithms) of the estimated counts are given in log_frequencies.txt
// Certain flags may be taken: 
// -r (resolution), size of the bins
// -t, max number of iterations
// -T, max time
// -S, number of samples at each bin to track, will be uniform over the configurations encountered
//          at that value

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <sstream>
#include <vector>
#include <iterator>
#include <random>
#include <vector>
#include <iomanip>
#include <map>

//*********************** GLOBAL VARIABLES **************************************************************
std::vector<std::vector<int>> neighbors;
std::vector<int> current_partition;
std::vector<int> histogram;
std::vector<double> log_frequency;
double current_conductance;
double ent_step;
int n, m, q;
int max_iters;
double bin_size;
int num_bins;
int MIN_MAX;
double FLAT_THRESH;

//*********************** FUNCTION DECLARATIONS **************************************************************

// W-L functions
void reset_histogram();
void rescale_frequencies();
int cond_to_ind(double cond);

// Conductance Calcuator
int cut_size(std::vector<int> partition);
int size(std::vector<int> partition);
double conductance(std::vector<int> partition);

// Helpers
void print_vector (std::vector<int> vec);
void print_vector (std::vector<double> vec);
void print(int val);
void print(double val);
double rng_double();

//*********************** MAIN PROGRAM ***********************************************************************

int main(int argc, char *argv[]) {
    // Input parameters
    max_iters = 10000000;
    bin_size = 0.01;
    MIN_MAX = 100;
    FLAT_THRESH = 0.95;

    // Read in.txt edgelist
    std::ifstream myfile ("in.txt"); // Filename
    if (!myfile.is_open()){
        std::cout << "Can't read in.txt \n";
        exit(0);
    }
    myfile >> n >>  m >> q; 
    // Instantiations
    for(int i = 0; i < n; i++){
      neighbors.push_back(std::vector<int>());
    }
    for (int i = 0; i < m; i++) {
        int e1, e2;
        myfile >> e1 >> e2;
        neighbors[e1].push_back(e2);
        neighbors[e2].push_back(e1);
    }
    // Find the range of conductances
    int max_degree = 0;
    for(int i = 0; i < n; i++){
      int degree = neighbors[i].size();
      max_degree = std::max(degree,max_degree);
    }
    num_bins = std::ceil(max_degree/bin_size) + 1;
    for(int i = 0; i < num_bins; i++){
      histogram.push_back(0);
      log_frequency.push_back(0);
    }
    // Generate initial partition (1/2 chance each)
    for(int i = 0; i < n;i++){
      if(rng_double() > 0.5){
        current_partition.push_back(0);
      }else{
        current_partition.push_back(1);
      }
    }
    // print_vector(current_partition);
    current_conductance = conductance(current_partition);
  // Main loop
  ent_step = 1;
  for(int t = 0; t < max_iters; t++){
    // Propose swap
    int swap_ind = std::floor(rng_double()*n);
    current_partition[swap_ind] = (current_partition[swap_ind]+1)%2; // Flip
    double proposed_conductance = conductance(current_partition);
    // print(cond_to_ind(proposed_conductance));
    double probAccept = exp(log_frequency[cond_to_ind(current_conductance)]- log_frequency[cond_to_ind(proposed_conductance)]);
    if(rng_double() < probAccept){
      // Keep change
      current_conductance = proposed_conductance;
    }else{
      // Reverse change
      current_partition[swap_ind] = (current_partition[swap_ind]+1)%2; // Flip
    }
    // Increment Histogram
    histogram[cond_to_ind(current_conductance)]++;
    log_frequency[cond_to_ind(current_conductance)]+=ent_step;
    // Check for flatness
    if(t % 100 == 0){
      int maxVal = 0;
      int minNonZero = 1000000; //inf
      int numNonZero = 0;
      for(int i = 0; i < num_bins; i++){
        maxVal = std::max(maxVal,histogram[i]);
        if(histogram[i] > 0){
          numNonZero++;
          minNonZero = std::min(minNonZero,histogram[i]);
        }
      }
        // cout << "Step: " << step << " Iterations: " << t << " Ratio: " << double(minNonZero)/maxVal << " # Found: " << numNonZero;
        // cout << " smallestCut: " << smallestCutSize << "/" << init_cut_size;
        // cout << " largestCut: " << largestCutSize << endl;
      if(maxVal > MIN_MAX && minNonZero > maxVal*FLAT_THRESH){
        std::cout << "Iterations: " << t << " ent_step: " <<ent_step << std::endl;
        ent_step /= 2;
        rescale_frequencies();
        reset_histogram();
        print(maxVal);
        // Write to file
        std::ofstream outfile("out.txt");
        outfile << bin_size << std::endl;
        for(int i = 0; i < num_bins; i++){
          outfile << log_frequency[i] << " ";
        }
        outfile.close();
      }
    }
  }
  rescale_frequencies();
  print_vector(log_frequency);
}

// W-L functions
void reset_histogram(){
  for(int i = 0; i < num_bins;i++){
    histogram[i] = 0;
  }
}
// add/subtract from the frequencies so that they actually reflect the log of the absolute number, using that
// they should sum to 2^(n-1).
void rescale_frequencies(){
  double max_log = *std::max_element(log_frequency.begin(), log_frequency.end());
  double exp_sum = 0;
  for(int i = 0; i < num_bins; i++){
    exp_sum += std::exp(log_frequency[i] - max_log);
  }
  for(int i = 0; i < num_bins; i++){
    log_frequency[i] = log_frequency[i] - max_log + (n-1)*std::log(2) - std::log(exp_sum);
  }
}

// Evalutations
int cut_size(std::vector<int> partition){
  int cut_number = 0;
  for(int i = 0; i < n; i++){
    for(int j:neighbors[i]){
      if(partition[i] != partition[j]){
        cut_number++;
      }
    }
  }
  return cut_number/2; // Undirected
}
int size(std::vector<int> partition){
  int total = std::accumulate(partition.begin(), partition.end(), 0);
  return std::min(n-total,total)*4; // Pick smaller size
}
double conductance(std::vector<int> partition){
  if(size(partition) == 0){
    return 0;
  }
  return cut_size(partition)/float(size(partition));
}
int cond_to_ind(double cond){
  return std::floor(cond/bin_size);
}

// Helper functions
double rng_double(){
  return static_cast<double>(std::rand()) / RAND_MAX;
}
void print_vector (std::vector<int> vec) {
  for (int n=0; n<vec.size(); ++n)
    std::cout << vec[n] << ' ';
  std::cout << '\n';
}
void print_vector (std::vector<double> vec) {
  for (int n=0; n<vec.size(); ++n)
    std::cout << vec[n] << ' ';
  std::cout << '\n';
}
void print(int val){
  std::cout<<val<<std::endl;
}
void print(double val){
  std::cout<<val<<std::endl;
}
