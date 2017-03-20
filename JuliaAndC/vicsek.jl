# gcc -Wall -shared -o calc.so -lm -fopenmp -fPIC calc.c
# export OMP_NUM_THREADS=4
# julia vicsek.jl
using ProgressMeter
unirand(a,b) = a + (b-a) * rand();

output_dir = "data"
if ! ispath(output_dir); mkdir(output_dir); end

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
        D = zeros(num_particle * num_particle);
        ccall((:ppdist, "./calc.so"), Void,( Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32, Float64), D, x, y, Int32(num_particle), L);

        # calculate mean angle
        mean_ang = zeros(num_particle);
        ccall((:cal_mean_angle, "./calc.so"), Void,( Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32, Float64), D, mean_ang, theta, Int32(num_particle), r);
      
        # update angle with noise
        theta = mean_ang .+ 2*pi*eta .* (rand(num_particle) .- 0.5);

        # update positions and velocities;
        vx = vc .* cos(theta); vy = vc .* sin(theta);
        x += vx * dt; y += vy * dt;

        # periodic boundary condition
        x = mod(x,L); y = mod(y,L);

        next!(progress);
    end

end;

#srand(1)
@time vicsek();
