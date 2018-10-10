public class Island {

	Individual[] population;

	public Island(Model model, int islandIdentity){

		// Constructor
		for (int i = 0; i < model.N; i++) {
			
			population[i] = new Individual(model, this);
			
		}

	}

	public Individual[] rankPopulation() {

		int N = population.length;
		Individual[] ranked_population = new Individual[];
		for (int i = 0; i < N; i++) {
			ranked_population[i] = new Individual(model, this);
		}
		double [] ranked_fitnesses = new double[N];
		double[] dummy = new double[N];
		for(int i = 0; i < N; i++) {

			dummy[i] = ranked_fitnesses[i] = population[i].fitness;
		}

		Arrays.sort(ranked_fitnesses);
		for(int i = 0; i < N; i++) {

			ranked_population[i] = population[find(dummy, ranked_fitnesses[i])];
		}
	}



}