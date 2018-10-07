import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.io.*;
import java.util.*;

public class player51 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
  private int evaluations_limit_;


	public player51()
	{
		// Initiate the pseudorandom number generator
		rnd_ = new Random();
	}


	public void setSeed(long seed)
	{
		// Set seed of algortihm's random process
		rnd_.setSeed(seed);
	}


	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;
		// Get evaluation properties
		Properties props = evaluation.getProperties();
    // Get evaluation limit
    evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
    boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
    boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
    boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm
    if(isMultimodal){
        // Do sth
    }else{
        // Do sth else
    }
  }


	public static void printArray(double matrix[][])
	{
		for (double[] row : matrix)
				System.out.println(Arrays.toString(row));
	}


	public void run()
	{
		int evals = 0;


		int N = 100;  // Number of individuals
		int D = 10;  // Dimension
		double F = 0.2;  // Scaling factor
		double Cr = 0.5;  // Crossover probability
		double [][] individuals = new double[N][D+1];  // 2D matrix used to keep track of individuals
		double [] champion = new double[D+1];

		//double [] initial = {3.656000004787539, 2.5495999982057698, -1.5296000398259961, 1.4695999840870468, 1.396000011255213, -1.908000003268257, 3.5015999931433224, -2.350400012612358, -0.38400001020493046, -2.0359999920987946};

		// Initial for Katsuur (score: 7.07) (9.987)
		//double [] initial = {4.710820431448943, 4.017982732015528, 4.903162145172683, 4.924815264346772, 2.6087568395233602, 4.509901451099045, 2.822437486928696, 4.635152856631775, 4.383806271805781, 4.770372449523396};
		double [] initial = {4.715029721575221, 4.016511185127753, 4.904475634824179, 4.923627614184904, 2.5962170543841725, 4.505780026656369, 2.824918466682541, 4.629462070803707, 4.388512044006422, 4.771607466228325};



		// Initialize matrix randomly and initial fitness evaluation
		for(int i=0; i<N; i++){
			for(int j=0; j<D; j++){
				//individuals[i][j] = rnd_.nextDouble()*5;
				individuals[i][j] = initial[j] + rnd_.nextDouble()*0.1;
				// if(rnd_.nextDouble() < 0.5){
				// 	individuals[i][j] = -1*individuals[i][j];
				// }
			}
			individuals[i][D] = (double) evaluation_.evaluate(Arrays.copyOfRange(individuals[i], 0, D));  // Assign some random "fitness" values
			evals++;
		}


		// Find first champion
		int index_best = 0;
		for(int i=0; i<N; i++){
			if(individuals[i][D] > individuals[index_best][D]){
				index_best = i;
			}
		}
		for(int i=0; i<D; i++){
			champion[i] = individuals[index_best][i];
		}
		champion[D] = individuals[index_best][D];


		while(evals < evaluations_limit_) {

			double [][] mutant_population = new double[N][D+1];
			double [][] trial_population = new double[N][D+1];

			for(int i=0; i<N; i++){
				double [] base_vector = individuals[rnd_.nextInt(N)];
				double [] parent1 = individuals[rnd_.nextInt(N)];
				double [] parent2 = individuals[rnd_.nextInt(N)];

				double [] perturbation_vector = new double[D+1];
				for(int j=0; j<D; j++){
					perturbation_vector[j] = F * (parent1[j] - parent2[j]);
					mutant_population[i][j] = base_vector[j] + perturbation_vector[j];
				}

				for(int j=0; j<D; j++){
					if(rnd_.nextDouble() < Cr){
						trial_population[i][j] = mutant_population[i][j];
					} else {
						trial_population[i][j] = individuals[i][j];
					}
				}

				int index = rnd_.nextInt(D);
				trial_population[i][index] = parent1[index];
			}

			for(int i=0; i<N; i++){
				trial_population[i][D] = (double) evaluation_.evaluate(Arrays.copyOfRange(trial_population[i], 0, D));
				evals++;
				if(trial_population[i][D] > individuals[i][D]){
					individuals[i] = trial_population[i];
				}
			}

			// Find champion
			index_best = 0;
			for(int i=0; i<N; i++){
				if(individuals[i][D] > individuals[index_best][D]){
					index_best = i;
				}
			}
			for(int i=0; i<D; i++){
				champion[i] = individuals[index_best][i];
			}
			champion[D] = individuals[index_best][D];

			if(evals >= evaluations_limit_-1){
				System.out.println(Arrays.toString(champion));
			}

		}  // End of while loop

	}  // End of run() function

}  // End of class
