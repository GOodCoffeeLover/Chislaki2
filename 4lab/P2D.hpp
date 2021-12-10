#pragma once
#include "Progonka.hpp"

#include <iomanip>
#include <cmath>
#include <iostream>
#include <functional>
#include <exception>
#include <vector>

using f_RxRxRtoR = std::function<double(double, double, double)>; // function f : RxRxR -> R
using f_RxRtoR = std::function<double(double, double)>; // function f : RxR -> R
using f_RtoR = std::function<double(double)>; // function f : R -> R


class P2D{
  public:
   
    P2D(){}
    ~P2D(){}
    void solve_ADM();
    void solve_FSM();
   
    void get_ans(std::vector<std::vector<std::vector<double>>>& ans_user);
    
    void print_ans_last_layer();

    void set_equation(double A, f_RxRxRtoR F);

    void set_left_border_condition (double a, double b, f_RxRtoR f);
    void set_right_border_condition(double a, double b, f_RxRtoR f);
    
    void set_upper_border_condition (double a, double b, f_RxRtoR f);
    void set_lower_border_condition(double a, double b, f_RxRtoR f);

    void set_zero_layer(f_RxRtoR f);

    void set_x_area(double l, double r, double step);
    void set_y_area(double l, double r, double step);
    void set_time(double T, double step);

    double MSE(f_RxRxRtoR);

  private:
    std::vector<std::vector<std::vector<double>> > ans{};

    // equation
    //   a*d^2u/dx^2 + a*d^2u/dx^2 = du/dt 
    struct equation{
      double a;
      f_RxRxRtoR f = [](double x, double y, double t) -> double{ return 0.0;}; 
    } equ;

    //border condition like
    // a*du(x0,y,t)/dx + b * u(x0,y,t) = f(x0,y, t) 
    struct border_condition{
      double a=0,b=1;
      f_RxRtoR f = [](double ax, double t)->double{ return ax*t;};
    } l_cond, r_cond, // along ox 
      u_cond, d_cond; // along oy

    // area of our equation for x-axis and time
    // from l to r with step = step
    struct area{
      double l, r, step;
    } x_axis, y_axis, time;

    f_RxRtoR zero_layer = [](double x, double y){ return x*y;};

};
