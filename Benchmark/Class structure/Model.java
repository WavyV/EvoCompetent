public class Model {

	// All global parameters
	public Island[] islands;
	public double mutationProbability;
	public int populationSize;
	public double sigma_init;

	public Model(double mutationProbability, int populationSize, double sigma_init,) {

		this.mutationProbability = mutationProbability;
		this.populationSize = populationSize;
		this.sigma_init = sigma_init;

		Island[0] = new PSO_island();	

	}

	
}