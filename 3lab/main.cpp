#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>

#include "EllipticEquation.hpp"


using namespace std;
struct arguments{
  bool error = false, time = false;
  double hx = 0.05, hy = 0.05;
};


void check_args(int argc, char** argv, arguments& args){
  
  for(int i=1; i<argc; ++i)
    if( string(argv[i]) == string("error") ){
      args.error = true;
      
    }else if( string(argv[i]) ==  string("time") ){
      args.time = true;

    }else if( (string(argv[i]) == string("hx")) && (i+1 < argc) ){
      char* p;
      double buf;
      buf = strtod(argv[i+1], &p);
      if(!*p){
        args.hx = buf;
      }else{
        cerr<<"Wrong Argumet"<<endl;
        cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step]"<<endl;
        exit(-1);
      }  
      i+=1;  
    
    }else if( (string(argv[i]) == string("hy")) && (i+1 < argc) ){
      char* p;
      double buf;
      buf = strtod(argv[i+1], &p);
      if(!*p){
        args.hy = buf;
      }else{
        cerr<<"Wrong Argumet"<<endl;
        cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step]"<<endl;
        exit(-1);
      }
      i+=1;    

    }else{
      cerr<<"Wrong Argumet"<<endl;
      cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step]"<<endl;
      exit(-1);
    }
  
}


int main(int argc, char* argv[]){
  
  arguments args;
  check_args(argc, argv, args);


  EllipticEquation solver;

  solver.set_equation([](double x, double y)->double{return 0.0;});

  solver.set_left_border_condition(0, 1, [](double y)->double{return y;});
  solver.set_right_border_condition(0, 1, [](double y)->double{return 1.0+y;});

  solver.set_upper_border_condition(0, 1, [](double x)->double{return 1.0+x;});
  solver.set_lower_border_condition(0, 1, [](double x)->double{return x;});


  solver.set_x_area(0, 1, args.hx);
  solver.set_y_area(0, 1, args.hy);

  
  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  //solvers: Libman Zeldel UpperRelax
  solver.solve(Zeldel, 0.8, 0.0000001, 10000);

  chrono::steady_clock::time_point end   = chrono::steady_clock::now();

  if(args.time)
    cerr<< "Time difference = " << chrono::duration_cast<chrono::microseconds>(end - begin).count()/1000.0 << "[ms]" << endl;
  
  if(!args.error)
    solver.print_ans();
  else{
    f_RxRtoR true_ans = [](double x, double y){ return x+y;};
    cerr<<"MSE = "<<solver.MSE(true_ans)<<endl;  
  }

  return 0;
}