#include <iostream>
#include <vector>
#include "iters.hpp"

using namespace std;



int main(){
  vector<double> v1={0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, v2(10, 1);
  cout<<v1<<endl;
  cout<<v2<<endl;


  v1 +=v2;
  cout<<prod(v1, v2)<<endl;
  
  // vector<vector<double>> matrix = 
  //   {{10, 1, 1}, 
  //    {2, 10, 1},
  //    {2, 2, 10}};
  // vector<double> b = {12, 13, 14}, ans = {};
  
  vector<vector<double>> matrix = 
    {{22, -3, -8, 7}, 
     {-8, -22, -4, -8},
     {8, -2, -18, 2},
     {7, 2, -9, -24}};
  vector<double> b = {-158, 254, -108, -24}, ans = {};
  
  libman(matrix, b, ans, 1.0, 0.00001);
  cout<<ans<<endl;

  return 0;
}