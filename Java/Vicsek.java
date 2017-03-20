import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;

public class Vicsek {
	// pairwise distance
	public static double dist(double x1, double x2, double y1, double y2, double L){
		double dx, dy;
		dx = Math.abs(x1 - x2);
		if (dx > L-dx)
			dx = L-dx;

		dy = Math.abs(y1 - y2);
		if (dy > L-dy)
			dy = L-dy;
		return Math.sqrt(dx*dx + dy*dy);
	}

	// uniform random number (a,b)
	public static double unirand(double a, double b){
		return a + (b-a) * Math.random();
	}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

	public static void main (String[] args) {
		int N=300, steps=500, near_count, count;
		double x[] = new double[N], y[] = new double[N], vx[] = new double[N], vy[] = new double[N], theta[] = new double[N];
		double D[][] = new double[N][N], mean_ang[] = new double[N], list_ang[] = new double[N];
		double L=5.0, vc=1.0, dt=0.1, eta=0.1, r=0.5, mean_sin, mean_cos;

		File newdir = new File("data");
		newdir.mkdir();

		// initial value
		for (int i=0; i<N; i++){
			theta[i] = unirand(-Math.PI, Math.PI);
			x[i] = unirand(0.0, L); //System.out.println(x[i]);
			y[i] = unirand(0.0, L);
			vx[i] = vc * Math.cos(theta[i]);
			vy[i] = vc * Math.sin(theta[i]);
		}

		// Start vicsek model
		for(int iter=0; iter<steps; iter++){
			// write x, y, vx, vy
			try{
				String filename = String.format("./data/%05d.csv", iter);
				File file = new File(filename);
				PrintWriter pw = new PrintWriter(file);
				for (int i=0; i<N; i++){
					pw.println(x[i] + ", " + y[i] + ", " + vx[i] + ", " + vy[i]);
				}
				pw.close();
			} catch (IOException ex) {
            ex.printStackTrace();
			}

			// pairwise distance
			for (int i=0; i<N-1; i++){
				for (int j=i+1; j<N; j++){
					D[i][j] = dist(x[i], x[j], y[i], y[j], L);
				}
			}
			for (int i=0; i<N-1; i++){
				for (int j=i+1; j<N; j++){
					D[j][i] = D[i][j];
				}
			}
			for (int i=0; i<N; i++){
				D[i][i] = 0.0;
			}


			//calculate mean angle
			for (int i=0; i<N; i++){
				near_count = 0;
				for (int j=0; j<N; j++){
					if (D[i][j] <= r){
						list_ang[near_count] = theta[j];
						near_count += 1;
					}
				}

				mean_sin = 0.0; mean_cos = 0.0;
				for (int j=0; j<near_count; j++){
					mean_sin += Math.sin(list_ang[j]);
					mean_cos += Math.cos(list_ang[j]);
				}
				mean_sin = (double) mean_sin / near_count;
				mean_cos = (double) mean_cos / near_count;

				mean_ang[i] = Math.atan2(mean_sin, mean_cos);
			}

			// update theta
			for (int i=0; i<N; i++){
				theta[i] = mean_ang[i] + 2*Math.PI*eta*unirand(-0.5, 0.5);
			}

			// update position and velocity
			for (int i=0; i<N; i++){
				vx[i] = vc * Math.cos(theta[i]); vy[i] = vc * Math.sin(theta[i]);
				x[i] += vx[i] * dt; y[i] += vy[i] * dt;
			}

			// periodic boundary condition
			for (int i=0; i<N; i++){
                x[i] = (x[i] + L) % L;
                y[i] = (y[i] + L) % L;
			}
		}
	}
}

