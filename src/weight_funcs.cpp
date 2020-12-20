#include <weight_funcs.h>

void dWpress(
	Eigen::Ref<const Eigen::Vector3d> d,
	double r, double H,
	Eigen::Vector3d& g
) {
	g = -45. / (M_PI * pow(H, 6)) * pow(H - r, 2) / r * d;
}


double Wpoly(double r, double H) {
	return 315. / (64. * M_PI * pow(H, 9)) * pow(H * H - r * r, 3);
}

void dWpoly(
	Eigen::Ref<const Eigen::Vector3d> d,
	double r, double H,
	Eigen::Vector3d& g
) {
	g = -945. / (32. * M_PI * pow(H, 9)) * pow(H * H - r * r, 2) * d;
}

double d2Wpoly(double r, double H) {
	return -945. / (32. * M_PI * pow(H, 9)) * (H * H - r * r) * (3. * H * H - 7. * r * r);
}

double d2Wvisco(double r, double H) {
	return 45. / (M_PI * pow(H, 6)) * (H - r);
}