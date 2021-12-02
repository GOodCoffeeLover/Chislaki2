#include "EllipticEquation.hpp"

#include "iters.hpp"
#include<iomanip>
using namespace std;



   
void EllipticEquation::get_ans(std::vector<std::vector<double>>& ans_user){
  ans_user = ans;
}


void EllipticEquation::set_equation(f_RxRtoR F){
  equ.f = F;
}

void EllipticEquation::set_left_border_condition (double a, double b, f_RtoR f){
  l_cond.a = a;
  l_cond.b = b;
  l_cond.f = f;

}

void EllipticEquation::set_right_border_condition(double a, double b, f_RtoR f){
  r_cond.a = a;
  r_cond.b = b;
  r_cond.f = f;
}

void EllipticEquation::set_upper_border_condition (double a, double b, f_RtoR f){
  u_cond.a = a;
  u_cond.b = b;
  u_cond.f = f;
}

void EllipticEquation::set_lower_border_condition(double a, double b, f_RtoR f){
  d_cond.a = a;
  d_cond.b = b;
  d_cond.f = f;
}


void EllipticEquation::set_x_area(double l, double r, double step){
  x_axis.l = l;
  x_axis.r = r;
  x_axis.step = step;
}
void EllipticEquation::set_y_area(double l, double r, double step){
  y_axis.l = l;
  y_axis.r = r;
  y_axis.step = step;
}


void EllipticEquation::solve(int solver, double omega=1.0, double eps = 0.0001, unsigned int iters=1000000){
  unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step +1,
               ny = (y_axis.r - y_axis.l)/y_axis.step +1;
  
  //cerr<<nx<<' '<<ny<<endl;

  vector<vector<double>> matrix(nx*ny, vector<double>(nx*ny, 0));
  vector<double> b(nx*ny, 0);
  //y=0
  //ny*j +i = uij

  //cout<<"lover"<<endl;
  for(unsigned int i=0; i<nx; i+=1){
    matrix[i][i]   += d_cond.b - d_cond.a/x_axis.step;
    matrix[i+1][i] += d_cond.a/x_axis.step;
    b[i] += d_cond.f(x_axis.l + double(i)*x_axis.step);
  }


  //cout<<"upper"<<endl;
  for(unsigned int i=0; i<nx; i+=1){
    matrix[nx*(ny-1)+i][nx*(ny-1)+i] = u_cond.b + u_cond.a/x_axis.step;
    matrix[nx*(ny-1)+i][nx*(ny-1)+i-1] = -u_cond.a/x_axis.step;
    b[nx*(ny-1)+i] = u_cond.f(x_axis.l + double(i)*x_axis.step);
  }


  //cout<<"left"<<endl;
  for(unsigned int j=1; j<ny-1; j+=1){
    matrix[nx*j+0][nx*j+0] = l_cond.b - l_cond.a/x_axis.step;
    matrix[nx*j+0][nx*j+1] = l_cond.a/x_axis.step;
    b[nx*j+0] = l_cond.f(y_axis.l + double(j)*y_axis.step);
  }


  //cout<<"right"<<endl;
  for(unsigned int j=1; j<ny-1; j+=1){
    matrix[nx*j+nx-1][nx*j+nx-1] = r_cond.b + r_cond.a/x_axis.step;
    matrix[nx*j+nx-1][nx*j+nx-2] = -r_cond.a/x_axis.step;
    b[nx*j+nx-1] = r_cond.f(y_axis.l + double(j)*y_axis.step);
  }
  
  //cout<<"inner"<<endl;
  for(unsigned int i=1; i<nx-1; i+=1){
    for(unsigned int j=1; j<ny-1; j+=1){
      matrix[nx*j+i][nx*j+i]     = -2.0/x_axis.step - 2.0/y_axis.step;
      matrix[nx*j+i][nx*j+i+1]   = 1.0/x_axis.step;
      matrix[nx*j+i][nx*j+i-1]   = 1.0/x_axis.step;
      matrix[nx*j+i][nx*(j+1)+i] = 1.0/y_axis.step;
      matrix[nx*j+i][nx*(j-1)+i] = 1.0/y_axis.step;
      b[nx*j+i] = equ.f(x_axis.l + double(i)*x_axis.step, y_axis.l + double(j)*y_axis.step);
    }
  }
  // for(int i =0; i<matrix.size(); ++i)
  //   cerr<<setw(2)<<i+1<<": "<<matrix[i]<<endl;
 
  //cout<<"solving"<<endl;
  
  vector<double>res{};


  switch(solver){
    case 1:
      libman(matrix, b, res, 1, eps, iters);
      break;
    
    case 2:
      zeldel(matrix, b, res, eps, iters);
      break;
    
    case 3:
      libman(matrix, b, res, omega, eps, iters);
      break;
    
  }


  ans.assign(nx+1, vector<double>(ny+1, 0.0));
  for(unsigned int i=0; i<=nx; ++i)
    for(unsigned int j=0; j<=ny; ++j)
      ans[i][j]=res[ny*j+i];
}

void EllipticEquation::print_ans(){
  unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step,
               ny = (y_axis.r - y_axis.l)/y_axis.step;
  
  for(unsigned int i=0; i<=nx; ++i)
    for(unsigned int j=0; j<=ny; ++j)
      cout<<x_axis.l + double(i)*x_axis.step<<' '
          <<y_axis.l + double(j)*y_axis.step<<' '
          <<ans[i][j]<<endl;
}

double EllipticEquation::MSE(f_RxRtoR true_ans){
  double error=0;
   unsigned int nx = (x_axis.r - x_axis.l)/x_axis.step,
               ny = (y_axis.r - y_axis.l)/y_axis.step;
  
  unsigned int size = ans.size()*ans[0].size();

  for(unsigned int i=0; i<=nx; ++i)
    for(unsigned int j=0; j<=ny; ++j){
      double ans_ij = true_ans(double(i)*x_axis.step, double(j)*y_axis.step);
      error += (ans_ij - ans[i][j])*(ans_ij - ans[i][j])/double(size);
    }
    return error;

}