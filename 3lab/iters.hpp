#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include<iostream>

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);

std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs);
std::vector<double> operator-(const std::vector<double>& lhs, const std::vector<double>& rhs);

std::vector<double> operator*(const double lhs, const std::vector<double>& rhs);
std::vector<double> operator*(const std::vector<double>& lhs, const double  rhs);
std::vector<double>& operator+=(std::vector<double>& lhs, const std::vector<double>& rhs);


double prod(const std::vector<double>& lhs, const std::vector<double>& rhs);
std::vector<double> abs(const std::vector<double>& v);


double norm(const std::vector<double>& l);
double norm(const std::vector<std::vector<double>>& matrix);


void libman(const std::vector<std::vector<double>>& matrix, std::vector<double> betha, std::vector<double>& ans, double omega, double eps, unsigned int);
void zeldel(const std::vector<std::vector<double>>& matrix, std::vector<double> betha, std::vector<double>& ans, double eps, unsigned int);