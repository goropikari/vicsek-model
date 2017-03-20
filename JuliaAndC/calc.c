#include <stdio.h>
#include <math.h>
#include <omp.h>

// prototype
double dist(double x1, double x2, double y1, double y2, double L);
void ppdist(double D[], double x[], double y[], int n, double L);
void cal_mean_angle(double D[], double mean_ang[], double theta[], int N, double r);
void cal_ith_mean_angle(double D[], double mean_ang[], double theta[], int N, double r, int i);

double dist(double x1, double x2, double y1, double y2, double L){
    double dx, dy;
    dx = fabs(x1 - x2);
    if (dx > L-dx)
        dx = L-dx;
    dy = fabs(y1 - y2);
    if (dy > L-dy)
        dy = L-dy;
    return sqrt(dx*dx + dy*dy);
}

void ppdist(double D[], double x[], double y[], int n, double L){
   int count;
   # pragma omp parallel for private(count)
   for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
         count = i*n + j;
         D[count] = dist(x[i], x[j], y[i], y[j], L);
      }
   }
}


void cal_mean_angle(double D[], double mean_ang[], double theta[], int N, double r){
	# pragma omp parallel for
    for (int i=0; i<N; i++){
        cal_ith_mean_angle(D, mean_ang, theta, N, r, i);
    }
}

void cal_ith_mean_angle(double D[], double mean_ang[], double theta[], int N, double r, int i){
	int near_count = 0;
    double list_ang[N];
	for (int j=0; j<N; j++){
		if (D[i*N+j] <= r){
			list_ang[near_count] = theta[j];
			near_count += 1;
		}
	}
	double mean_sin=0.0, mean_cos=0.0;
	for (int j=0; j<near_count; j++){
		mean_sin += sin(list_ang[j]);
		mean_cos += cos(list_ang[j]);
	}
   
	mean_sin = mean_sin / near_count;
	mean_cos = mean_cos / near_count;
   
	mean_ang[i] = atan2(mean_sin, mean_cos);
}
