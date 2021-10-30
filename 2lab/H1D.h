#ifndef __H1D_LIB__
#define __H1D_LIB__

#include <functional>
#include <cmath>
#include <iostream>
#include <exception>
#include <vector>
using f_RxRtoR = std::function<double(double, double)>; // function f : RxR -> R
using f_RtoR = std::function<double(double)>; // function f : R -> R

#define FIRST_APROX_SECOND_COND 0
#define SECOND_APROX_SECOND_COND 1

#define EXPLICIT_SCHEME 0
#define IMPLICIT_SCHEME 1


class H1D{
  public:
   
    H1D(){}
    ~H1D(){}
    void solve(int APROX_SECOND_COND, int SHEME_TYPE);
    void get_ans(std::vector<std::vector<double>>& ans_user);
    
    void set_equation(double a, double b, double c, f_RxRtoR F);

    void set_left_border_condition (double a, double b, f_RtoR f);
    void set_right_border_condition(double a, double b, f_RtoR f);

    void set_start_layer_cond(f_RtoR F0, f_RtoR F1);

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
    struct border_condition{
      double a=0,b=1;
      f_RtoR f = [](double t)->double{ return 0.0;};
    } l_cond, r_cond;

    //condition for zero layer
    f_RtoR F0 = [](double x)->double{return std::sin(M_PI * 2.0 * x);};

    f_RtoR F1 = [](double x)->double{return std::sin(M_PI * 2.0 * x);};

    // area of our equation for x-axis and time
    // from l to r with step = step
    struct area{
      double l, r, step;
    } x_axis, time;

};

#endif