#include "P1D.h"
#include "Progonka.hpp"

using namespace std;

void P1D::set_equation(double a, double b, double c, f_RxRtoR F){
  if(a==0.0)
    throw logic_error("a param of equ is zero\n it's not p1d equ");
    
  ans.clear();
  equ.a=a;
  equ.b=b;
  equ.c=c;
  equ.f=F;
}

void P1D::set_left_edge_cond (double a, double b, f_RtoR f){
  if(a+b==0.0)
    throw logic_error("Left condition is uqeation 0 + 0 = f(t)");
    ans.clear();
    l_cond.a=a;
    l_cond.b=b;
    l_cond.f = f;
}
void P1D::set_right_edge_cond(double a, double b, f_RtoR f){
  if(a+b==0.0)
    throw logic_error("Right condition is uqeation 0 + 0 = f(t)");
    ans.clear();
    r_cond.a=a;
    r_cond.b=b;
    r_cond.f = f;
}

void P1D::set_zero_layer_cond(f_RtoR n_F){
  ans.clear();
  P1D::F=n_F;
}

void P1D::set_x_area(double l, double r, double step){
  ans.clear();
  x_axis.l = l;
  x_axis.r = r;
  x_axis.step = step; // \h
}
void P1D::set_t_area(double l, double r, double step){
  ans.clear();
  time.l = l;
  time.r = r;
  time.step = step; // \tau
}

void P1D::solve(double teta){ 
// teta == 0 : explicit (yavnaya)
// teta == 1 : implicit (neyavnaya)
// teta == 1/2 : Krank-Nikolson
  ans.clear();
  
  ans.push_back({});
  for(double x=x_axis.l; x<=x_axis.r; x+=x_axis.step)
    ans[0].push_back(F(x));

  for(double t=time.l+time.step; t<=time.r+0.0000000001; t+=time.step){
    
    vector<vector<double>> matrix{};
    
    if(l_cond.a==0)
      matrix.push_back({0,l_cond.b,0,l_cond.f(t)});
    else
      matrix.push_back(
        {0,
         2*pow(equ.a, 2)/x_axis.step + x_axis.step/time.step - equ.c*x_axis.step 
         - l_cond.b/l_cond.a*(2*pow(equ.a, 2) - equ.b*x_axis.step),
         -2*pow(equ.a, 2)*x_axis.step,
         x_axis.step/time.step*ans.back()[0] - l_cond.f(t)*(2*pow(equ.a, 2) - equ.b*x_axis.step)/l_cond.a});

    
    // j var for iteration in last layer
    for(int j = 1; j < (x_axis.r - x_axis.l)/x_axis.step; ++j ){
      matrix.push_back({
        -teta*pow(equ.a,2.0)/pow(x_axis.step,2.0) + equ.b/(2.0*x_axis.step), //j-1
        1.0/time.step + 2.0*teta*pow(equ.a,2.0)/pow(x_axis.step,2) - equ.c,  //j
        -teta*pow(equ.a,2.0)/pow(x_axis.step,2.0) - equ.b/(2.0*x_axis.step), //j+1
        ans.back()[j]/time.step + (1.0-teta)*pow(equ.a,2.0)/pow(x_axis.step,2.0)
        *(ans.back()[j+1] - 2*ans.back()[j] + ans.back()[j-1]) + equ.f(x_axis.l+j*x_axis.step, t) //d
      });


    }

    
    if(r_cond.a==0)
      matrix.push_back({0,r_cond.b,0,r_cond.f(t)});
    else
      matrix.push_back(
        {-2*pow(equ.a, 2)*x_axis.step,
         2*pow(equ.a, 2)/x_axis.step + x_axis.step/time.step - equ.c*x_axis.step 
         + r_cond.b/r_cond.a*(2*pow(equ.a, 2) + equ.b*x_axis.step),
         0,
         x_axis.step/time.step*ans.back().back() - r_cond.f(t)*(2*pow(equ.a, 2) + equ.b*x_axis.step)/r_cond.a});
    ans.push_back({});
    Progonka(matrix, ans.back());

  }

  

}


void P1D::get_ans(vector<vector<double>>& ans_user){
  ans_user = ans;
}

