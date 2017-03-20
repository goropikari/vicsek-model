from numpy.random import *
import numpy as np
import csv, os
def dist(x1,x2,L):
	dx = abs(x1[0] - x2[0]) % L
	dx = min([dx, L-dx])
	dy = abs(x1[1] - x2[1]) % L
	dy = min([dy, L-dy])
	return np.sqrt(dx**2 + dy**2)

def ppdist(x ,y, num_particle, L):
    D = np.zeros([num_particle, num_particle])
    for i in range(num_particle-1):
        for j in range(i+1, num_particle):
			D[i][j] = dist([ x[i], y[i] ], [ x[j], y[j] ], L);
    return D


def cal_mean_angle(D, theta, num_particle, r):
    meanangle = np.zeros(num_particle)
    for iter2 in range(num_particle):
        condition = D[iter2,:] < r;
        near_particles_angles = theta[condition];
        meanangle[iter2] = \
        np.arctan2( \
    	    np.mean( np.sin(near_particles_angles) ), \
    	    np.mean( np.cos(near_particles_angles) ) \
        );
    return meanangle;
	
output_dir = 'data'
if not os.path.exists(output_dir): os.makedirs(output_dir)

# parameter
# L:             system size
# num_particle: the number of particles
# vc:             constant speed
# dt:             time step
# steps:         the number of run steps
# eta:             Order
# r:             interaction radius
L = 5.0; num_particle = 300; vc = 1.0;
dt = 0.1; steps = 500; eta = 0.1; r = 0.5;


# initial condition
x = uniform(0.0, L, num_particle); y = uniform(0.0, L, num_particle); theta = uniform(0.0, L, num_particle)
vx = vc * np.cos(theta); vy = vc * np.sin(theta);


# start fish swimming!
#progress = Progress(steps);
for iter in range(steps):
	# output configuration
	filename = "./data/{0:05d}.csv".format(iter)
	data = np.transpose([x, y, vx, vy])
	f = open(filename, 'w')
	writer = csv.writer(f)
	writer.writerows(data)
	f.close()
	
	# pairwise distance
	# D[i,j] is a distance between i th and j th particles.
	D = ppdist(x, y, num_particle, L);
	D += np.transpose(D);

	# calculate mean angle
	meanangle = np.zeros(num_particle);
	meanangle = cal_mean_angle(D, theta, num_particle, r);
	
	# update angle with noise
	theta = meanangle + 2 * np.pi * eta * uniform(-0.5, 0.5, num_particle);

	# update positions and velocities;
	vx = vc * np.cos(theta); vy = vc * np.sin(theta);
	x += vx * dt; y += vy * dt;

	# periodic boundary condition
	x = x % L; y = y % L;
