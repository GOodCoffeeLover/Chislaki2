#include "P2D.hpp"
#include "Progonka.hpp"




   
void P2D::get_ans(std::vector<std::vector<std::vector<double>>>& ans_user){
  ans_user = ans;
}


void P2D::set_equation(double A, f_RxRxRtoR F){
  equ.f = F;
}

void P2D::set_left_border_condition (double a, double b, f_RxRtoR f){
  l_cond.a = a;
  l_cond.b = b;
  l_cond.f = f;
}

void P2D::set_right_border_condition(double a, double b, f_RxRtoR f){
  r_cond.a = a;
  r_cond.b = b;
  r_cond.f = f;
}

void P2D::set_upper_border_condition (double a, double b, f_RxRtoR f){
  u_cond.a = a;
  u_cond.b = b;
  u_cond.f = f;
}

void P2D::set_lower_border_condition(double a, double b, f_RxRtoR f){
  d_cond.a = a;
  d_cond.b = b;
  d_cond.f = f;
}

void P2D::set_zero_layer(f_RxRtoR F){
  zero_layer = F;
}

void P2D::set_x_area(double l, double r, double step){
  x_axis.l = l;
  x_axis.r = r;
  x_axis.step = step;
}

void P2D::set_y_area(double l, double r, double step){
  y_axis.l = l;
  y_axis.r = r;
  y_axis.step = step;
}

void P2D::set_time(double T, double step){
  time.l = 0;
  time.r = T;
  time.step = step;
}


// Alternating Direction Method
void P2D::solve_ADM(){
  unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step + 1, 
               ny = (y_axis.r - y_axis.l)/y_axis.step + 1, 
               nt = (  time.r -   time.l)/time.step   + 1;
  ans.assign(nt, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
  double x,y,t;
  for(unsigned int j=0; j<ny; ++j)
    for(unsigned int i=0; i<nx; ++i){
      x=x_axis.l+double(i)*x_axis.step;
      y=y_axis.l+double(j)*y_axis.step;

      ans[0][i][j]=zero_layer(x,y);
    }

  // for k = 0..nt-1
  //   for j = 0..ny
  //    buf[j][i] = solve u(xi, yj, tk+1/2):
  //   for i = 0..nx   
  //    ans[k+1][i][j] = solve u(buff[j][i], yj, tk+1)
  //

  std::vector<std::vector<double>> buff(ny, std::vector<double>(nx, 0.0));
  for(unsigned int k=0; k<nt-1; ++k){
    t = time.l + double(k)*time.step;
    // std::cerr<<t<<std::endl;

    for(unsigned int j=1; j<ny-1; ++j){
      y=y_axis.l+double(j)*y_axis.step;
      std::vector<std::vector<double>> mtrx{};
      
      mtrx.push_back({0, l_cond.b, 0, l_cond.f(y,t+time.step/2.0)});

      for(unsigned int i=1; i<nx-1; ++i)
        mtrx.push_back({
          -equ.a/x_axis.step*x_axis.step,
          2.0/time.step + 2.0*equ.a/x_axis.step*x_axis.step,
          -equ.a/x_axis.step*x_axis.step,

          2.0/time.step*ans[k][i][j] + equ.a / y_axis.step * y_axis.step 
          * (ans[k][i][j+1] -2.0*ans[k][i][j] + ans[k][i][j-1])});
      
      mtrx.push_back({0, r_cond.b, 0, r_cond.f(y,t+time.step/2.0)});  

      Progonka(mtrx, buff[j]);

    }

    // for(unsigned i=0; i<ny; ++i){
    //   x=x_axis.l+double(i)*x_axis.step;
    //   buff[0][i] = d_cond.f(x,t + time.step/2.0);
    //   buff[ny-1][i] = u_cond.f(x,t + time.step/2.0);
    // }

    for(unsigned int i=1; i<nx-1; ++i){
      x=x_axis.l+double(i)*x_axis.step;
      std::vector<std::vector<double>> mtrx{};
      
      mtrx.push_back({0.0, d_cond.b, 0.0, d_cond.f(x,t + time.step)});

      for(unsigned int j=1; j<ny-1; ++j)
        mtrx.push_back({
          -equ.a/y_axis.step*y_axis.step,
          
          2.0/time.step + 2.0*equ.a/y_axis.step*y_axis.step,
          
          -equ.a/y_axis.step*y_axis.step,

          2.0/time.step * buff[j][i] + equ.a / x_axis.step * x_axis.step 
          * (buff[j][i+1] -2.0 * buff[j][i] + buff[j][i-1])});
      
      mtrx.push_back({0.0, u_cond.b, 0.0, u_cond.f(x,t + time.step)});  

      Progonka(mtrx, ans[k+1][i]);
    }
    
    for(unsigned j=0; j<ny; ++j){
      y=y_axis.l+double(j)*y_axis.step;
      ans[k+1][0][j] = l_cond.f(y,t+time.step);
      ans[k+1][nx-1][j] =r_cond.f(y,t+time.step);
    }
  
  }
}

//Fractional Steps Method
void P2D::solve_FSM(){
   unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step + 1, 
                ny = (y_axis.r - y_axis.l)/y_axis.step + 1, 
                nt = (  time.r -   time.l)/time.step   + 1;
  ans.assign(nt, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
  double x,y,t;
  for(unsigned int j=0; j<ny; ++j)
    for(unsigned int i=0; i<nx; ++i){
      x=x_axis.l+double(i)*x_axis.step;
      y=y_axis.l+double(j)*y_axis.step;

      ans[0][i][j]=zero_layer(x,y);
    }

  // for k = 0..nt-1
  //   for j = 0..ny
  //    buf[j][i] = solve u(xi, yj, tk+1/2):
  //   for i = 0..nx   
  //    ans[k+1][i][j] = solve u(buff[j][i], yj, tk+1)
  //

  std::vector<std::vector<double>> buff(ny, std::vector<double>(nx, 0.0));
  for(unsigned int k=0; k<nt-1; ++k){
    t = time.l + double(k)*time.step;
    // std::cerr<<t<<std::endl;

    for(unsigned int j=1; j<ny-1; ++j){
      y=y_axis.l+double(j)*y_axis.step;
      std::vector<std::vector<double>> mtrx{};
      
      mtrx.push_back({0, l_cond.b, 0, l_cond.f(y,t+time.step/2.0)});

      for(unsigned int i=1; i<nx-1; ++i)
        mtrx.push_back({
          -equ.a/x_axis.step*x_axis.step,
          2.0/time.step + 2.0*equ.a/x_axis.step*x_axis.step,
          -equ.a/x_axis.step*x_axis.step,

          2.0/time.step*ans[k][i][j]});
      
      mtrx.push_back({0, r_cond.b, 0, r_cond.f(y,t+time.step/2.0)});  

      Progonka(mtrx, buff[j]);

    }

    // for(unsigned i=0; i<ny; ++i){
    //   x=x_axis.l+double(i)*x_axis.step;
    //   buff[0][i] = d_cond.f(x,t + time.step/2.0);
    //   buff[ny-1][i] = u_cond.f(x,t + time.step/2.0);
    // }

    for(unsigned int i=1; i<nx-1; ++i){
      x=x_axis.l+double(i)*x_axis.step;
      std::vector<std::vector<double>> mtrx{};
      
      mtrx.push_back({0.0, d_cond.b, 0.0, d_cond.f(x,t + time.step)});

      for(unsigned int j=1; j<ny-1; ++j)
        mtrx.push_back({
          -equ.a/y_axis.step*y_axis.step,
          
          2.0/time.step + 2.0*equ.a/y_axis.step*y_axis.step,
          
          -equ.a/y_axis.step*y_axis.step,

          2.0/time.step * buff[j][i]});
      
      mtrx.push_back({0.0, u_cond.b, 0.0, u_cond.f(x,t + time.step)});  

      Progonka(mtrx, ans[k+1][i]);
    }
    
    for(unsigned j=0; j<ny; ++j){
      y=y_axis.l+double(j)*y_axis.step;
      ans[k+1][0][j] = l_cond.f(y,t+time.step);
      ans[k+1][nx-1][j] =r_cond.f(y,t+time.step);
    }
  
  }
}

void P2D::print_ans_last_layer(){
  unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step,
               ny = (y_axis.r - y_axis.l)/y_axis.step;
  
  for(unsigned int i=0; i<=nx; ++i)
    for(unsigned int j=0; j<=ny; ++j)
      std::cout<<x_axis.l + double(i)*x_axis.step<<' '
          <<y_axis.l + double(j)*y_axis.step<<' '
          <<ans.back()[i][j]
          <<std::endl;
}

double P2D::MSE(f_RxRxRtoR true_ans){
  double error=0;
  unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step + 1, 
               ny = (y_axis.r - y_axis.l)/y_axis.step + 1, 
               nt = (  time.r -   time.l)/time.step   + 1;
  
  double x,y,t;

  for(unsigned int k=0; k<nt; ++k) 
    for(unsigned int j=0; j<ny; ++j) 
      for(unsigned int i=0; i<nx; ++i){
        x=x_axis.l+double(i)*x_axis.step;
        y=y_axis.l+double(j)*y_axis.step;
        t=time.l+double(k)*time.step;
        error += (ans[k][i][j] - true_ans(x,y,t))*(ans[k][i][j] - true_ans(x,y,t));
      }
  error /= double(ans.size()*ans[0].size()*ans[0][0].size());
 
  return error;

}