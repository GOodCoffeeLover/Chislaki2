#include <iostream>
#include <cmath>
#include <chrono>
#include "P2D.hpp"

using namespace std;


struct arguments{
  bool error = false, time = false, FSM = true;
  double hx = M_PI/16.0, hy = M_PI/16.0, ht = 0.001;
};


void check_args(int argc, char** argv, arguments& args){
  
  for(int i=1; i<argc; ++i)
    if( string(argv[i]) == string("error") ){
      args.error = true;
      
    }else if( string(argv[i]) ==  string("time") ){
      args.time = true;

    }else if( string(argv[i]) ==  string("FSM") ){
      args.FSM = true;

    }else if( string(argv[i]) ==  string("ADM") ){
      args.FSM = false;

    }else if( (string(argv[i]) == string("hx")) && (i+1 < argc) ){
      char* p;
      double buf;
      buf = strtod(argv[i+1], &p);
      if(!*p){
        args.hx = buf;
      }else{
        cerr<<"Wrong Argumet"<<endl;
        cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step] [ht <t_step>] [ADM] [FSM]"<<endl;
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
        cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step] [ht <t_step>] [ADM] [FSM]"<<endl;
        exit(-1);
      }
      i+=1;    

    }else if( (string(argv[i]) == string("ht")) && (i+1 < argc) ){
      char* p;
      double buf;
      buf = strtod(argv[i+1], &p);
      if(!*p){
        args.ht = buf;
      }else{
        cerr<<"Wrong Argumet"<<endl;
        cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step] [ht <t_step>] [ADM] [FSM]"<<endl;
        exit(-1);
      }
      i+=1;    

    }else{
      cerr<<"Wrong Argumet"<<endl;
      cerr<<"Usage: "<<argv[0]<<" [error] [time] [hx <x_step>] [hy y_step] [ht <t_step>] [ADM] [FSM]"<<endl;
      exit(-1);
    }
  
}
int main(int argc, char* argv[]){
  P2D solver;

  arguments args;
  check_args(argc, argv, args);

  double a=1, mu1=1, mu2=1, T=0.1;

  solver.set_equation(a, [](double x, double y, double t) -> double{ return 0.0;});

  solver.set_left_border_condition(0, 1, 
    [a, mu1, mu2](double y, double t)->double{ 
      return cos(mu2*y) * exp( -(mu1*mu1 + mu2*mu2) * a * t);
    });
  solver.set_right_border_condition(0, 1,  
    [a, mu1, mu2](double y, double t)->double{return (mu1==1? -1: 1)*cos(mu2*y)*exp(-(mu1*mu1 + mu2*mu2)*a*t);});

  solver.set_lower_border_condition(0, 1, 
    [a, mu1, mu2](double x, double t)->double{return cos(mu1*x)*exp(-(mu1*mu1 + mu2*mu2)*a*t);});
  solver.set_upper_border_condition(0, 1,  
    [a, mu1, mu2](double x, double t)->double{return (mu2==1? -1: 1)*cos(mu1*x)*exp(-(mu1*mu1 + mu2*mu2)*a*t);});

  solver.set_zero_layer([mu1, mu2](double x, double y)->double{return cos(mu1*x)*cos(mu2*y);});

  solver.set_x_area(0, M_PI, args.hx);
  solver.set_y_area(0, M_PI, args.hy);
  solver.set_time(T, args.ht);

  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  if(args.FSM){
    solver.solve_FSM();
  }else{
    solver.solve_ADM();
  }

  chrono::steady_clock::time_point end   = chrono::steady_clock::now();

  if(args.time)
    cerr<< "Time difference = " << chrono::duration_cast<chrono::microseconds>(end - begin).count()/1000.0 << "[ms]" << endl;


  if(args.error)
    cerr<<"Error = "<<solver.MSE([a, mu1, mu2](double x, double y, double t){ return cos(mu1*x)*cos(mu2*y)*exp(-(mu1*mu1 + mu2*mu2)*a*t);})<<endl;
  else
    solver.print_ans_last_layer();

  return 0;

}