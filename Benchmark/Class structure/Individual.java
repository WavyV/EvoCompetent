public class Individual {


	public double[] values;
	public double[] sigma;
	public double[][] alpha;
	public double fitness;
	public Island homeland;
	public Model model;

	public Individual(Model model, Island homeland) { // args must be of order int N (pop size), boolean Cor, Boolean Uncor
		

		rnd_ = new Random();
		this.homeland = homeland;
		this.model = model;

		for (int i = 0; i < N; i++) {

			values[i] = rnd_.newGaussian() * model.sigma_init;
		}
	}

	public void setFitness(newFitness) {

		fitness = newFitness;

	}

	public void mutate_Gaussian {

		if Math.random() < 

	}
}