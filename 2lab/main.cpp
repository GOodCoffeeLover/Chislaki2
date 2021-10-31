#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <exception>

#include "H1D.h"

using namespace std;

const int steps=20;

// 0     1     2     3 4    5     6
// a.out xstep tstep a apox sheme error
int main(int argc, char ** argv){

  try{
    H1D equ;
  
    double a=(argc>=4?stod(argv[3]):0.5),b=0,c=0;
    f_RxRtoR f_eq = [](double x,double t)->double{return 0.0;}; 
    equ.set_equation(a,b,c,f_eq);

    double alpha=0, beta=1;
    f_RtoR f_l_cond = [a](double t)->double{return -sin(a*t);};
    equ.set_left_border_condition(alpha, beta, f_l_cond);
    
    double gamma=0, delta=1;
    f_RtoR f_r_cond = [a](double t)->double{return sin(a*t);};
    equ.set_right_border_condition(gamma, delta, f_r_cond);

    f_RtoR f_zero_layer = [](double x)->double{return sin(x);};
    f_RtoR f_diff = [a](double x)->double{return -a*cos(x);};
    
    equ.set_start_layer_cond(f_zero_layer, f_diff);

    double xl=0.0, xr = M_PI, xstep=(argc>=2?stod(argv[1]):M_PI/steps);
    equ.set_x_area(xl, xr, xstep);

    double tl=0.0, tr = 1.0, tstep=(argc>=3?stod(argv[2]):1.0/steps);
    equ.set_t_area(tl, tr, tstep);

    equ.solve(SECOND_APROX_SECOND_COND, IMPLICIT_SCHEME);
    equ.solve(
      (argc>=5?stoi(argv[4]):SECOND_APROX_SECOND_COND), 
      (argc>=6?stoi(argv[5]):IMPLICIT_SCHEME) 
    );
    vector<vector<double>> ans{};

    equ.get_ans(ans);

    if(argc>=7 && string(argv[6])=="error"){
      f_RxRtoR real_ans_f=[a](double x, double t){ return sin(x-a*t);}; 

      double error=0.0;
      for(int j=0; j<ans.back().size(); ++j){
        error+=abs(real_ans_f(xl+xstep*j, tr) - ans.back()[j]);
      }
      error/=ans.back().size();
      cout<<"MAE = "<<error<<endl;
    }else{
      for(int k=0; k<ans.size(); ++k)
        for(int j=0; j<ans[k].size(); ++j)
          cout<<(tl+tstep*k)<<' '<<(xl+xstep*j)<<' '<<ans[k][j]<<endl;
    }

  }catch (exception& e){
    cout << e.what() <<endl;
    return 1;
  }
  


  return 0;
}