#include "H1D.h"
#include "Progonka.hpp"



using namespace std;

void H1D::set_equation(double a, double b, double c, f_RxRtoR F){
  if(a==0.0)
    throw logic_error("a param of equ is zero\n it's not H1D equ");
  if(a<0.0)
    throw logic_error("a param of equ is less then zero");
    
  ans.clear();
  equ.a=a;
  equ.b=b;
  equ.c=c;
  equ.f=F;
}

void H1D::set_left_border_condition (double a, double b, f_RtoR f){
  if((abs(a)+abs(b))==0.0)
    throw logic_error("Left condition is uqeation 0 + 0 = f(t)");
    ans.clear();
    l_cond.a=a;
    l_cond.b=b;
    l_cond.f = f;
}
void H1D::set_right_border_condition(double a, double b, f_RtoR f){
  if((abs(a)+abs(b))==0.0)
    throw logic_error("Right condition is uqeation 0 + 0 = f(t)");
    ans.clear();
    r_cond.a=a;
    r_cond.b=b;
    r_cond.f = f;
}

void H1D::set_start_layer_cond(f_RtoR n_F0, f_RtoR n_F1 ){
  ans.clear();
  H1D::F0=n_F0;
  H1D::F1=n_F1;
}

void H1D::set_x_area(double l, double r, double step){
  ans.clear();
  x_axis.l = l;
  x_axis.r = r;
  x_axis.step = step; // \h
}
void H1D::set_t_area(double l, double r, double step){
  ans.clear();
  time.l = l;
  time.r = r;
  time.step = step; // \tau
}

void H1D::solve(int APROX_SECOND_COND, int SHEME_TYPE){ 
  if(pow(equ.a, 2.0)*time.step/pow(x_axis.step, 2.0) > 0.5 && SHEME_TYPE==EXPLICIT_SCHEME)
    throw logic_error("a^2*\\tau/h^2 > 1/2 for explicit method");

  ans.clear();
  double xj, tk;
  
  ans.push_back({});
  for(int j = 0; j <= (x_axis.r - x_axis.l)/x_axis.step; ++j ){
    xj=x_axis.l + x_axis.step*j; 
    ans[0].push_back(F0(xj));
  }

  ans.push_back({});
  switch(APROX_SECOND_COND){

    case FIRST_APROX_SECOND_COND:
      for(int j = 0; j <= (x_axis.r - x_axis.l)/x_axis.step; ++j ){
        xj=x_axis.l + x_axis.step*j; 
        ans[1].push_back(
          F0(xj) + 
          F1(xj) * time.step);
      }
      break;

    case SECOND_APROX_SECOND_COND:
      int F0xx;
      for(int j = 0; j <= (x_axis.r - x_axis.l)/x_axis.step; ++j ){
        xj=x_axis.l + x_axis.step*j; 
        F0xx = (F0(xj+x_axis.step) - 2.0*F0(xj) + F0(xj-x_axis.step))/pow(x_axis.step, 2.0);
        ans[1].push_back( 
          F0(xj) + 
          F1(xj) * time.step + 
          pow(equ.a,2.0) * F0xx * pow(time.step, 2.0) / 2.0);
      }
      
      break;


    default:
      ans.clear();
      throw logic_error("unexpected value of APROX_SECOND_COND");
  }



  for(int k=2; k<(time.r -time.l)/time.step; k+=1){
    tk=time.l + k*time.step;
    
    vector<vector<double>> matrix{};
    
    if(l_cond.a==0)
      matrix.push_back({0,l_cond.b,0,l_cond.f(tk)});
    else
      matrix.push_back(
        {0,
         2*pow(equ.a, 2)/x_axis.step + x_axis.step/time.step - equ.c*x_axis.step 
         - l_cond.b/l_cond.a*(2*pow(equ.a, 2) - equ.b*x_axis.step),
         -2*pow(equ.a, 2)*x_axis.step,
         x_axis.step/time.step*ans[k-1][0] - l_cond.f(tk)*(2*pow(equ.a, 2) - equ.b*x_axis.step)/l_cond.a});

    
    // j var for iteration in last layer
    switch(SHEME_TYPE){
      
      case EXPLICIT_SCHEME: 
        for(int j = 1; j < (x_axis.r - x_axis.l)/x_axis.step; ++j ){
          xj=x_axis.l + x_axis.step*j; 
          matrix.push_back({
            equ.b/(2.0*x_axis.step), //j-1
            1.0/ pow(time.step, 2.0) + 1.0/time.step - equ.c,  //j
            - equ.b/(2.0*x_axis.step), //j+1
            
            ans[k-1][j]/time.step + 
            ans[k-1][j-1]*pow(equ.a/x_axis.step, 2.0) + 
            ans[k-1][j]*(2.0/pow(time.step, 2.0) -2.0*pow(equ.a/x_axis.step, 2.0) ) + 
            ans[k-1][j+1]*pow(equ.a/x_axis.step, 2.0) + 
            -ans[k-2][j]/pow(time.step, 2.0) +
            equ.f(xj, tk) //d
          });
        }

        break;
      
      case IMPLICIT_SCHEME: 

        for(int j = 1; j < (x_axis.r - x_axis.l)/x_axis.step; ++j ){
          xj=x_axis.l + x_axis.step*j; 
          
          matrix.push_back({
            - pow(equ.a/x_axis.step, 2.0) + equ.b/(2.0*x_axis.step), //j-1
            2.0*pow(equ.a/x_axis.step, 2.0) + 1.0/pow(time.step, 2.0) + 1.0/time.step - equ.c,  //j
            - pow(equ.a/x_axis.step, 2.0) - equ.b/(2.0*x_axis.step), //j+1
            ans[k-1][j]/time.step + 
            ans[k-1][j]*2.0/pow(time.step, 2.0) -
            ans[k-2][j]/pow(time.step, 2.0) + 
            equ.f(x_axis.l+j*x_axis.step, tk) //d
          });
        }
        
        break;

      default:
        ans.clear();
        throw logic_error("unexpected value of SHEME_TYPE");

    }

    
    if(r_cond.a==0)
      matrix.push_back({0,r_cond.b,0,r_cond.f(tk)});
    else
      matrix.push_back(
        {-2*pow(equ.a, 2)*x_axis.step,
         2*pow(equ.a, 2)/x_axis.step + x_axis.step/time.step - equ.c*x_axis.step 
         + r_cond.b/r_cond.a*(2*pow(equ.a, 2) + equ.b*x_axis.step),
         0,
         x_axis.step/time.step*ans[k-1].back() - r_cond.f(tk)*(2*pow(equ.a, 2) + equ.b*x_axis.step)/r_cond.a});
    ans.push_back({});
    Progonka(matrix, ans.back());

  }

  

}


void H1D::get_ans(vector<vector<double>>& ans_user){
  ans_user = ans;
}

