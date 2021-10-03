#ifndef __P1D_LIB__
#define __P1D_LIB__

#include <functional>
#include <cmath>
#include <iostream>
#include <exception>
#include <vector>
using f_RxRtoR = std::function<double(double, double)>; // function f : RxR -> R
using f_RtoR = std::function<double(double)>; // function f : R -> R


class P1D{
  public:
   
    P1D(){}
    ~P1D(){}
    void solve(double teta);
    void get_ans(std::vector<std::vector<double>>& ans_user);
    
    void set_equation(double a, double b, double c, f_RxRtoR F);

    void set_left_edge_cond (double a, double b, f_RtoR f);
    void set_right_edge_cond(double a, double b, f_RtoR f);

    void set_zero_layer_cond(f_RtoR F);

    void set_x_area(double l, double r, double step);
    void set_t_area(double l, double r, double step);


  private:
    std::vector<std::vector<double> > ans{};

    // equation
    // du/dt = a^2 * d^2u/dx^2 + b * du/dx + c * u + f(x,t) 
    struct equation{
      double a=1, b=0, c=0;
      f_RxRtoR f = [](double x, double t) -> double{ return 0.0;}; 
    } equ;

    //edge condition
    // a*du(x0,t)/dx + b * u(x0,t) = f(t) 
    struct edge_cond{
      double a=0,b=1;
      f_RtoR f = [](double t)->double{ return 0.0;};
    } l_cond, r_cond;

    //condition for zero layer
    f_RtoR F = [](double x)->double{return std::sin(M_PI * 2.0 * x);};

    // area of our equation for x-axis and time
    // from l to r with step = step
    struct area{
      double l, r, step;
    } x_axis, time;

};

#endif