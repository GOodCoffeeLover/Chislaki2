#ifndef __EE_LIB__
#define __EE_LIB__

#include <functional>
#include <cmath>
#include <iostream>
#include <exception>
#include <vector>
#include "iters.hpp"
using f_RxRtoR = std::function<double(double, double)>; // function f : RxR -> R
using f_RtoR = std::function<double(double)>; // function f : R -> R

#define Libman 1
#define Zeldel 2
#define UpperRelax 3

class EllipticEquation{
  public:
   
    EllipticEquation(){}
    ~EllipticEquation(){}
    void solve(int, double, double, unsigned int);
   
    void get_ans(std::vector<std::vector<double>>& ans_user);
    
    void print_ans();

    void set_equation(f_RxRtoR F);

    void set_left_border_condition (double a, double b, f_RtoR f);
    void set_right_border_condition(double a, double b, f_RtoR f);
    
    void set_upper_border_condition (double a, double b, f_RtoR f);
    void set_lower_border_condition(double a, double b, f_RtoR f);


    void set_x_area(double l, double r, double step);
    void set_y_area(double l, double r, double step);

    double MSE(f_RxRtoR);

  private:
    std::vector< std::vector<double> > ans{};

    // equation
    //   d^2u/dx^2 + d^2u/dx^2 = f(x,y) 
    struct equation{
      f_RxRtoR f = [](double x, double y) -> double{ return 0.0;}; 
    } equ;

    //border condition like
    // a*du(x0,y)/dx + b * u(x0,y) = f(x0,y) 
    struct border_condition{
      double a=0,b=1;
      f_RtoR f = [](double t)->double{ return t;};
    } l_cond, r_cond, u_cond, d_cond;

    // area of our equation for x-axis and time
    // from l to r with step = step
    struct area{
      double l, r, step;
    } x_axis, y_axis;

};

#endif