#include "iters.hpp"
#include<iomanip>

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v){
    for(double d: v)
      os<<std::setw(3)<<d<<' ';
    
    return os;
}

double prod(const std::vector<double>& lhs, const std::vector<double>& rhs){
  if(lhs.size() != rhs.size())
    std::cerr<<"Warning: Iputed vectors different lenghts: "<<lhs.size()<<" and "<<rhs.size()<<std::endl;
  unsigned int n = std::min(lhs.size(), rhs.size());
  double ans=0.0;
  for(unsigned int i=0; i<n; ++i)
    ans +=lhs[i]*rhs[i];

  return ans;
}

std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs){
  
  if(lhs.size() != rhs.size())
    std::cerr<<"Warning: Iputed vectors different lenghts: "<<lhs.size()<<" and "<<rhs.size()<<std::endl;
  
  unsigned int n = std::min(lhs.size(), rhs.size());
  std::vector<double> res(n);
  
  for(unsigned int i=0; i<n; ++i)
    res[i] =lhs[i]+rhs[i];

  return res;
}


std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs){
  
  if(lhs.size() != rhs.size())
    std::cerr<<"Warning: Iputed vectors different lenghts: "<<lhs.size()<<" and "<<rhs.size()<<std::endl;
  
  unsigned int n = std::min(lhs.size(), rhs.size());
  std::vector<double> res(n);
  
  for(unsigned int i=0; i<n; ++i)
    res[i] =lhs[i]-rhs[i];

  return res;
}

std::vector<double> operator*(const double lhs, const std::vector<double>& rhs){
  std::vector<double> res(rhs);
  for(unsigned int i=0; i<rhs.size(); ++i)
    res[i] =lhs*rhs[i];

  return res;
}

std::vector<double> operator*(const std::vector<double>& lhs, const double  rhs){
  
  return rhs*lhs;
}

std::vector<double>& operator+=(std::vector<double>& lhs, const std::vector<double>& rhs){
  if(lhs.size() != rhs.size())
    std::cerr<<"Warning: Iputed vectors different lenghts: "<<lhs.size()<<" and "<<rhs.size()<<std::endl;
  unsigned int n = std::min(lhs.size(), rhs.size());
  for(unsigned int i=0; i<n; ++i)
    lhs[i] +=rhs[i];

  return lhs;
}


std::vector<double> abs(const std::vector<double>& v){
  std::vector<double> vv(v);
  for(double& d: vv)
    d=std::abs(d);
  return vv;
}

double norm(const std::vector<double>& l){
  std::vector<double> ones(l.size(), 1);
  return prod(abs(l), ones);
}

double norm(const std::vector<std::vector<double>>& matrix){
  std::vector<double> ones(matrix.size(), 1);
  double m = norm(matrix[0]);

  for(unsigned int i=1; i<matrix.size(); ++i){
    m = std::max(m, norm(matrix[i]));
  }

  return m;
}


void libman(const std::vector<std::vector<double>>& matrix, std::vector<double> betha, std::vector<double>& ans, double omega, double eps, unsigned int iters){

  for(unsigned int i = 0; i<matrix.size(); ++i)
    if(matrix.size() != matrix[i].size())
      throw std::logic_error("non-quadric matrix");

  unsigned int n = matrix.size();
  std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.0));
  std::vector<double> buff(n, 0.0);

  ans.assign(n, 0.0);
  
  for(unsigned int i = 0; i<n; ++i)
    for(unsigned int j = 0; j<n; ++j)
      if(i==j){
        betha[i] /=matrix[i][i];
        ans[i] =betha[i];
      }else 
        alpha[i][j]=-matrix[i][j]/matrix[i][i];
  

  double a=norm(alpha);
  if(a>=1)
    std::cerr<<"Warning: ||alpha|| >= 1 \nAlpha: "<<a<<std::endl;
  unsigned int iter;
  for(iter = 0; iter<iters; ++iter){
    for(unsigned int i = 0; i<n; ++i){
      buff[i] = prod(alpha[i], ans) + betha[i];
    }
    
    if((a!= 1.0 && std::pow(a, iter+1)/(1.0 - a)*norm(ans-buff) <= eps) || (norm(ans-buff) <= eps)){
      std::cerr<<"Iter finish:"<<iter<<std::endl;
      break;
    }
    ans=omega*buff + (1-omega)*ans;
  }
}


void zeldel(const std::vector<std::vector<double>>& matrix, std::vector<double> betha, std::vector<double>& ans, double eps, unsigned int iters){

  for(unsigned int i = 0; i<matrix.size(); ++i)
    if(matrix.size() != matrix[i].size())
      throw std::logic_error("non-quadric matrix");

  unsigned int n = matrix.size();
  std::vector<std::vector<double>> alpha(n, std::vector<double>(n, 0.0));

  ans.assign(n, 0.0);
  
  for(unsigned int i = 0; i<n; ++i)
    for(unsigned int j = 0; j<n; ++j)
      if(i==j){
        betha[i] /=matrix[i][i];
        ans[i] =betha[i];
      }else 
        alpha[i][j]=-matrix[i][j]/matrix[i][i];
  

  std::vector<double> prev(ans);
  double a=norm(alpha);
  if(a>=1)
    std::cerr<<"Warning: ||alpha|| >= 1 \nAlpha: "<<a<<std::endl;
  for(unsigned int iter = 0; iter<iters; ++iter){
    for(unsigned int i = 0; i<n; ++i){
      ans[i] = prod(alpha[i], ans) + betha[i];
    }

    // if(iter%10==0)
    //   std::cerr<<iter<<" : "<<norm(ans-prev)<<std::endl;
    if((a!= 1.0 && std::pow(a, iter+1)/(1.0 - a)*norm(ans-prev) <= eps) || (norm(ans-prev) <= eps)){
      std::cerr<<"Iter finish:"<<iter<<std::endl;
      break;
    }
    prev=ans;
  }

}