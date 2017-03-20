using PyPlot, ProgressMeter;

output_dir = "pic"
if ! ispath(output_dir); mkdir(output_dir); end

function pic()
    L = 5.0 # system size
    steps = 500;
    progress = Progress(steps);
    for i in 0:steps-1
        filename = @sprintf("./data/%05d.csv", i);
        data = readcsv(filename);
        x = data[:,1]; y = data[:,2]; vx = data[:,3]; vy = data[:,4];
        plot(x,y,"b.");
        quiver(x,y,vx,vy);
        axis("square"); axis([0, L, 0, L]);
        picname = @sprintf("./pic/%05d.png", i);
        savefig(picname, dpi=60);

        clf();
        next!(progress);
    end
end
pic()
gc();
