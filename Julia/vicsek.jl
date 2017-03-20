using ProgressMeter # If you didn't install ProgressMeter, do Pkg.add("ProgressMeter")
unirand(a,b) = a + (b-a) * rand(); # define fn returning uniform random number in (a,b)

# make directory to store the positions and the velocities data.
output_dir = "data"
if ! ispath(output_dir); mkdir(output_dir); end

# distane between i th particle and j th particle at periodic boundary condition
function dist(x1::Array{Float64,1} ,x2::Array{Float64,1} ,L::Float64);
    # x[1]: x position of a particle
    # x[2]: y position of a particle
    # L: system size
    dx = abs(x1[1] - x2[1]) % L;
    dx = minimum([dx; L-dx]);
    dy = abs(x1[2] - x2[2]) % L;
    dy = minimum([dy; L-dy]);
    return sqrt(dx^2 + dy^2)
end

# pairwise distance
function ppdist(x::Array{Float64,1} ,y::Array{Float64,1}, num_particle::Int64, L::Float64);
    D = zeros(num_particle, num_particle);
    for i in 1:num_particle-1
            for j in i+1:num_particle
                    D[i,j] = dist([ x[i]; y[i] ], [ x[j]; y[j] ], L);
            end
    end

    return D;
end


# calculate mean angle
function cal_mean_angle(D::Array{Float64,2}, theta::Array{Float64,1}, num_particle::Int64, r::Float64)
    meanangle = zeros(num_particle)
    for iter2 in 1:num_particle
        truth = D[:,iter2] .< r;
        near_particles_angles = theta[truth];
        meanangle[iter2] =
        atan2(
            mean( sin(near_particles_angles) ),
            mean( cos(near_particles_angles) )
        );
    end
    return meanangle;
end

function vicsek()
    # L:             system size
    # num_particle: the number of particles
    # vc:             constant speed
    # dt:             time step
    # steps:         the number of run steps
    # eta:             Order
    # r:             interaction radius
    const L = 5.0; const num_particle = 300; const vc = 1.0;
    const dt = 0.1; const steps = 500; const eta = 0.1; const r = 0.5;

    # initial position of the particles
    # x_i and y_i corespond to the i th positon of particle.
    x = rand(num_particle) * L; y = rand(num_particle) * L;

    # initial angle (-pi, pi) of the particles;
    theta = [unirand(-pi, pi) for i in 1:num_particle];
    #theta = unirand(-pi, pi, num_particle);
    

    # initial velocity
    vx = vc .* cos(theta); vy = vc .* sin(theta);

    # start fish swimming!
    progress = Progress(steps);
    for iter in 0:steps-1;
        # output configuration
        filename = @sprintf("./data/%05d.csv",iter);
        writecsv(filename, [x y vx vy])
      
        # pairwise distance
        # D[i,j] is a distance between i th and j th particles.
        D = ppdist(x, y, num_particle, L);
        D += D';

        # calculate mean angle
        meanangle = zeros(num_particle);
        meanangle = cal_mean_angle(D, theta, num_particle, r);
      
        # update angle with noise
        theta = meanangle .+ 2*pi*eta .* (rand(num_particle) .- 0.5);

        # update positions and velocities;
        vx = vc .* cos(theta); vy = vc .* sin(theta);
        x += vx * dt; y += vy * dt;

        # periodic boundary condition
        #x = (x .+ L) % L; y = (y .+ L) % L;
        x = mod(x,L); y = mod(y,L);

        next!(progress);
    end  

end;
@time vicsek();
