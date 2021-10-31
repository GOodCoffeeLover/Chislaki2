#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <exception>

#include "P1D.h"



using namespace std;

int main(int argc, char ** argv){

  P1D equ;
  
  double a=(argc>=5?stod(argv[4]):0.1),b=0,c=0;
  f_RxRtoR f_eq = [](double x,double t)->double{return 0.0;}; 
  equ.set_equation(sqrt(a),b,c,f_eq);

  double alpha=0, beta=1;
  f_RtoR f_l_cond = [a](double t)->double{return exp(-a*t);};
  equ.set_left_edge_cond(alpha, beta, f_l_cond);
  
  double gamma=0, delta=1;
  f_RtoR f_r_cond = [a](double t)->double{return -exp(-a*t);};
  equ.set_right_edge_cond(gamma, delta, f_r_cond);

  f_RtoR f_zero_layer = [](double x)->double{return cos(x);};
  equ.set_zero_layer_cond(f_zero_layer);

  //cout<<"h = "<<stod(argv[1])<<endl;
  double xl=0.0, xr = M_PI, xstep=(argc>=2?stod(argv[1]):M_PI/4.0);
  equ.set_x_area(xl, xr, xstep);

  //cout<<"\\tau = "<<stod(argv[2])<<endl;
  double tl=0.0, tr = 1.0, tstep=(argc>=3?stod(argv[2]):0.3333);
  equ.set_t_area(tl, tr, tstep);

  try{
    double teta = (argc>=4 && (stod(argv[3])>=0.0 && stod(argv[3])<=1.0)?stod(argv[3]):0.5);
    equ.solve(teta);

  }catch (exception& e)
  {
    cout << e.what() <<endl;
    return 1;
  }

  vector<vector<double>> ans{};

  equ.get_ans(ans);

  if(argc>=6 && string(argv[5])=="error"){
    f_RxRtoR real_ans_f=[a](double x, double t){ return exp(-4.0*M_PI*M_PI*a*t)*sin(2.0*M_PI*x);}; 

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

  


  return 0;
}