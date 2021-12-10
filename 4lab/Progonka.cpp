#include "Progonka.hpp"
//matrix               ||   mtrix format
//b1 c1  0  0 | d1  ||   { 0, b1, c1, d1}
//a2 b2 c2  0 | d2  ||   {a2, b2, c2, d2}
// 0 a3 b3 c3 | d3  ||   {a3, b3, c3, d3}
// 0  0 a4 b4 | d4  ||   {a4, b4, 0, d4}

void Progonka(const std::vector<std::vector<double>> &mtrix, std::vector<double>& ans){
    ans.clear();
    std::vector<double> P(mtrix.size()+1, 0.0), 
                 Q(mtrix.size()+1, 0.0); 
                 ans.assign(mtrix.size(), 0.0);
  
  for(int i=1; i<P.size(); ++i)
    //Pi = -ci/(bi + ai * P(i-1))
    P[i]=-mtrix[i-1][2]/(mtrix[i-1][1] + mtrix[i-1][0]*P[i-1]); 
  
  for(int i=1; i<Q.size(); ++i)
    //Qi = (di - ai*Q(i-1)) / (bi + ai * P(i-1))
    Q[i] = (mtrix[i-1][3] - mtrix[i-1][0]*Q[i-1])/(mtrix[i-1][1] + mtrix[i-1][0]*P[i-1]);
  
  //xn = qn
  ans[ans.size()-1] = Q[Q.size()-1];
  for(int i=ans.size()-1; i>0; --i)
    //x(i) = Qi + Pi * x(i+1)
    ans[i-1]=Q[i] + P[i]*ans[i];
}