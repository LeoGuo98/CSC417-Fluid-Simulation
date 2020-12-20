#ifndef WEIGHT_FUNCS_H
#define WEIGHT_FUNCS_H

#include <Eigen/dense>

double Wpoly(double r, double H);

void dWpoly(
	Eigen::Ref<const Eigen::Vector3d> d,
	double r, double H,
	Eigen::Vector3d& g
);

double d2Wpoly(double r, double H);

void dWpress(
	Eigen::Ref<const Eigen::Vector3d> d,
	double r, double H,
	Eigen::Vector3d& g
);

double d2Wvisco(double r, double H);

#endif