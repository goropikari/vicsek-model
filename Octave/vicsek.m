clear;clc;
tic
# pairwise distane when periodic boundary condition
function f = dist(x1 ,x2 ,L);
    # x(1): x position of a particle
    # x(2): y position of a particle
    # L: system size
    dx = mod(abs(x1(1) - x2(1)), L);
    dx = min([dx; L-dx]);
    dy = mod(abs(x1(2) - x2(2)), L);
    dy = min([dy; L-dy]);
    f = sqrt(dx^2 + dy^2);
end

function g = ppdist(x, y, num_particle, L)
    D = zeros(num_particle, num_particle);
    for i = 1:num_particle-1
        for j = i+1:num_particle
            x1 = [x(i); y(i)];, x2 = [x(j); y(j)];
            D(i,j) = dist(x1, x2, L);
        end
    end
    g = D;
end


# calculate mean angle
function h = cal_mean_angle(D, theta, num_particle, r)
    meanangle = zeros(num_particle,1);
    for iter2 = 1:num_particle
        condition = D(:, iter2) < r;
        near_particles_angles = theta(condition);
        meanangle(iter2) = ...
        atan2( ...
            mean( sin(near_particles_angles) ), ...
            mean( cos(near_particles_angles) ) ...
        );
    end
    h = meanangle;
end

# L:            system size
# num_particle: the number of particles
# vc:           constant speed
# dt:           time step
# steps:        the number of run steps
# eta:          Order
# r:            interaction radius
L = 5.0; num_particle = 300; vc = 1.0;
dt = 0.1; steps = 500; eta = 0.1; r = 0.5;

# initial position of the particles
# x_i and y_i corespond to the i th positon of particle.
x = rand(num_particle,1) * L; y = rand(num_particle,1) * L;

# initial angle (-pi, pi) of the particles;
theta = -pi .+ 2*pi.*(rand(num_particle,1) .- 0.5);

# initial velocity
vx = vc .* cos(theta); vy = vc .* sin(theta);

# start fish swimming!
mkdir("data");
for iter = 0:steps-1;
    # output configuration
    filename = sprintf("./data/%05d.csv",iter);
    csvwrite(filename, [x y vx vy])
    
    # pairwise distance
    # D[i,j] is a distance between i th and j th particles.
    D = zeros(num_particle);
    D = ppdist(x, y, num_particle, L);
    D += D';

    # calculate mean angle
    meanangle = zeros(num_particle,1);
    meanangle = cal_mean_angle(D, theta, num_particle, r);
    
    # update angle with noise
    theta = meanangle .+ 2*pi*eta .* (rand(num_particle,1) .- 0.5);

    # update positions and velocities;
    vx = vc * cos(theta); vy = vc * sin(theta);
    x += vx * dt; y += vy * dt;

    # periodic boundary condition
    x = mod(x,L); y = mod(y,L);

end 
toc
