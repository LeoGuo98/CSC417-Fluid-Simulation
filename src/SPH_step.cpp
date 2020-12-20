#include <SPH_step.h>
#include <vector>
#include <weight_funcs.h>
#include <coefs.h>
#include <iostream>
#include <collision_response.h>

void SPH_step(Particles& particles, double dt, Coef& coef) {


	double
		density,
		d2Cs, // Color surface Laplacian
		r; // positions distance
	Eigen::Vector3d
		d,
		f_pressure, // pressure 
		f_visco, // viscosity 
		f_surface, // surface tension 
		f_g, // gravity 
		grad_press, // W_pressure gradient
		grad_poly, // W_poly gradient
		dCs; // Color surface normal 

	double min_x = particles.position.col(0).minCoeff();
	double min_y = particles.position.col(1).minCoeff();
	double min_z = particles.position.col(2).minCoeff();

	double max_x = particles.position.col(0).maxCoeff();
	double max_y = particles.position.col(1).maxCoeff();
	double max_z = particles.position.col(2).maxCoeff();

	int cell_x = (max_x - min_x) / coef.H + 1;
	int cell_y = (max_y - min_y) / coef.H + 1;
	int cell_z = (max_z - min_z) / coef.H + 1;

	std::vector<std::vector<std::vector<std::vector<int>>>> grid(cell_x, std::vector<std::vector<std::vector<int>>>(cell_y, std::vector<std::vector<int>>(cell_z, std::vector<int>())));
	int grid_x, grid_y, grid_z;

	for (size_t i = 0; i < particles.position.rows(); i++) {
		Eigen::Vector3d pos = particles.position.row(i);
		grid_x = (pos(0) - min_x) / coef.H;
		grid_y = (pos(1) - min_y) / coef.H;
		grid_z = (pos(2) - min_z) / coef.H;
		grid[grid_x][grid_y][grid_z].push_back(i);
	}



	for (grid_x = 0; grid_x < cell_x; grid_x++) {
		for (grid_y = 0; grid_y < cell_y; grid_y++) {
			for (grid_z = 0; grid_z < cell_z; grid_z++) {
				for (int i : grid[grid_x][grid_y][grid_z]) {
					density = 0;
					for (int count_x = -1; count_x < 2; count_x++) {

						if (grid_x + count_x < 0 || grid_x + count_x >= cell_x) continue;

						for (int count_y = -1; count_y < 2; count_y++) {

							if (grid_y + count_y < 0 || grid_y + count_y >= cell_y) continue;

							for (int count_z = -1; count_z < 2; count_z++) {

								if (grid_z + count_z < 0 || grid_z + count_z >= cell_z) continue;

								for (int j : grid[grid_x + count_x][grid_y + count_y][grid_z + count_z]) {
									d = particles.position.row(i) - particles.position.row(j);
									r = sqrt(d.dot(d));
									if (r > coef.H) continue;
									density += Wpoly(r, coef.H);
								}
							}
						}
					}

					particles.density(i) = density * coef.MASS;
					particles.pressure(i) = coef.STIFFNESS * (particles.density(i) - coef.RHO_IDEAL);
				}
			}
		}
	}


	for (grid_x = 0; grid_x < cell_x; grid_x++) {
		for (grid_y = 0; grid_y < cell_y; grid_y++) {
			for (grid_z = 0; grid_z < cell_z; grid_z++) {

				for (int i : grid[grid_x][grid_y][grid_z]) {

					// forces and surface tension

					f_pressure.setZero();
					f_surface.setZero();
					f_g << 0, coef.G* particles.density(i), 0;
					f_visco.setZero();
					dCs.setZero();
					d2Cs = 0;

					for (int count_x = -1; count_x < 2; count_x++) {

						if (grid_x + count_x < 0 || grid_x + count_x >= cell_x) continue;

						for (int count_y = -1; count_y < 2; count_y++) {

							if (grid_y + count_y < 0 || grid_y + count_y >= cell_y) continue;

							for (int count_z = -1; count_z < 2; count_z++) {

								if (grid_z + count_z < 0 || grid_z + count_z >= cell_z) continue;

								for (int j : grid[grid_x + count_x][grid_y + count_y][grid_z + count_z]) {
									d = particles.position.row(i) - particles.position.row(j);
									r = sqrt(d.dot(d));
									if (r > coef.H) continue;
									if (r > 0.) {
										dWpress(d, r, coef.H, grad_press);
										f_pressure -= coef.MASS * (particles.pressure(i) + particles.pressure(j)) / (2. * particles.density(j)) * grad_press;
										Eigen::Vector3d diff_v = particles.velocity.row(j) - particles.velocity.row(i);
										f_visco += d2Wvisco(r, coef.H) * diff_v / particles.density(j);
									}
									dWpoly(d, r, coef.H, grad_poly);
									dCs += particles.density(j) * grad_poly;
									d2Cs += d2Wpoly(r, coef.H) / particles.density(j);
								}
							}
						}
					}

					f_visco *= coef.VISCOSITY * coef.MASS;
					dCs *= coef.MASS;
					d2Cs *= coef.MASS;



					if (dCs.norm() > coef.SUFRACE_THRESH) {
						f_surface = -coef.SUFRACE_TENSION * d2Cs * dCs.normalized();
					}
					particles.acceleration.row(i) = (f_pressure + f_visco + f_surface + f_g) / particles.density(i);
					particles.velocity.row(i) += dt * particles.acceleration.row(i);
					particles.position.row(i) += dt * particles.velocity.row(i);
				}
			}
		}
	}
}
