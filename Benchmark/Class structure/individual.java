public class Individual {


	public double[] values;
	public double[] sigma;
	public double[][] alpha;
	public double fitness;
	private double[][] covMatrix;

	public Individual(Model ea, int N, double sigma_init, boolean Cor, boolean Uncor) { // args must be of order int N (pop size), boolean Cor, Boolean Uncor
		

		rnd_ = new Random();
		covMatrix = covMatrix;

		for (int i = 0; i < N; i++) {

			values[i] = rnd_.newGaussian() * sigma_init
		}

		if (Cor || Uncor) {

			sigma[i] = rnd_.newGaussian() * sigma_init
			if (Uncor) {

				for (int j = 0; j < N; j++) {

					alpha[i][j] = 0

				}
			}
		}		
	}

	public void setFitness(newFitness) {

		fitness = newFitness;

	}

	public void mutate_Gaussian {

		if Math.random() < 

	}
}