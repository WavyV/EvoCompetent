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

	// This function will return the index of some value in an array
	public static int find(double[] array, double value)
	{
		int index = 0;
		for(int i=0; i<array.length; i++){
				 if(array[i] == value){
						 index = i;
						 break;
				 }
		}
		return index;
	}

	// This function will rank the population (low index=low fitness, high index=high fitness)
	public static double[][] rank_population(double[][] population)
	{
		int N = population.length;
		int M = population[0].length;

		double [][] ranked_population = new double[N][M];
		double [] ranked_fitnesses = new double[N];
		double [] dummy = new double[N];

		for(int i=0; i<N; i++){
			ranked_fitnesses[i] = population[i][M-1];
			dummy[i] = population[i][M-1];
		}
		Arrays.sort(ranked_fitnesses);

		for(int i=0; i<N; i++){
			ranked_population[i] = population[find(dummy, ranked_fitnesses[i])];
		}

		return ranked_population;
	}

	// This function prints a matrix to the console
	public static void printArray(double matrix[][])
	{
		for (double[] row : matrix)
				System.out.println(Arrays.toString(row));
	}


	public void run()
	{

		int evals = 0;

    int N = 100;  // Number of individuals (make sure N/2 is even)
    int D = 10;  // Dimension of problem
    double [][] individuals = new double[N][D+2];  // 2D matrix used to keep track of individuals first 10 columns are the problem values, 11th is the fitness
    double sigma_init = 2.5;  // Sigma of Gaussian used in initialization
    double tau = 0.3;  // Sigma of Gaussian used in mutation
		double epsilon = 0.00001;  // Minimum value for sigma


    // Initialize matrix randomly and initial fitness evaluation
    for(int i=0; i<N; i++){
      for(int j=0; j<D+1; j++){
				individuals[i][j] = rnd_.nextGaussian()*sigma_init;
				if(j == D){
					individuals[i][j] = sigma_init*rnd_.nextGaussian();
				}
				// } else {
				// 	individuals[i][j] = initial[j] + rnd_.nextGaussian()*0.5;
				// }
      }
  		individuals[i][D+1] = (double) evaluation_.evaluate(Arrays.copyOfRange(individuals[i], 0, D));
			evals++;  // Important to increase evals after every evaluation, otherwise an error will occur
    }

    // Initial ranking of the population. Rank[0]=bad and Rank[N]=good.
    individuals = rank_population(individuals);

		// In this loop we go through the parent selection, crossover, mutation, survivor selection
		// until we reach the maximum number of allowed evaluations.
		while(evals < evaluations_limit_) {

      // Parent selection
			// Here we want to pick parents from the last half of the individuals matrix (that is the best half)
			// dummy_array keeps track of which parents we have already picked
			// parent1 and parent2 keep the information of the current parents we have picked
			// Once both parent1/2 have been picked crossover is done
			// Then both children are put in the children matrix
			// Process is repeated until entire best half of the population has mated
      int parent_number = 1;
			int [] dummy_array = new int[N];
      double [] parent1 = new double[D+2];
      double [] parent2 = new double[D+2];
      double [][] children = new double[N/2][D+2];

      for(int i=0; i<N/2; i++){
				// We want to pick N/2 parents in total
        int random_index = 0;
        boolean picked = false;
        do {
					// Keep picking random indices until we have found a parent that has not mated yet
          random_index = rnd_.nextInt(N/2) + N/2;
          if(dummy_array[random_index] == 0){
            dummy_array[random_index] = 1;
            picked = true;
          }
        } while(picked == false);

				// If we have picked 2 parents, then go do crossover, else pick another parent
        if(parent_number == 1){
          parent1 = individuals[random_index];
          parent_number++;
        } else {
          parent2 = individuals[random_index];
          parent_number = 1;

          // Create children using one-point crossover
          int random_split = rnd_.nextInt(D) + 1;
          for(int j=0; j<random_split; j++){
            children[i-1][j] = parent1[j];
            children[i][j] = parent2[j];
          }
          for(int j=random_split; j<(D+2); j++){
            children[i-1][j] = parent2[j];
            children[i][j] = parent1[j];
          }
        }
      }


			// Update sigmas
			for(int i=0; i<children.length; i++){
				children[i][D] = children[i][D] * Math.exp(tau*rnd_.nextGaussian());
				if(children[i][D] < epsilon){
					children[i][D] = epsilon;
				}
			}


      // Apply mutation, simple fixed sigma Gaussian mutation with 50% probability on each gene
      for(int i=0; i<children.length; i++){
        for(int j=0; j<D; j++){
          if(rnd_.nextDouble() < 1){
            children[i][j] += rnd_.nextGaussian()*children[i][D];
          }
        }
      }

      // Check fitness of children
      for(int i=0; i<children.length; i++){
        children[i][D+1] = (double) evaluation_.evaluate(Arrays.copyOfRange(children[i], 0, D));
				evals++;
      }

      // Select survivors
			// First put all the previous generation and the children together in one matrix
			// Then rank the matrix based on fitness
			// The select the best N to form the new generation
      double [][] total_population = new double[individuals.length+children.length][D+2];
      for(int i=0; i<individuals.length; i++){
        total_population[i] = individuals[i];
      }
      for(int i=0; i<children.length; i++){
        total_population[i+individuals.length] = children[i];
      }
      total_population = rank_population(total_population);

			// Selecting the new generation at the end of the array to have it ranked properly
      int index = individuals.length-1;
      for(int i=total_population.length-1; i >= children.length; i--){
        individuals[index] = total_population[i];
        index--;
      }


			if(evals >= evaluations_limit_-1){
				System.out.println(Arrays.toString(individuals[N-1]));
			}

    }  // End of EA loop

  }  // End of run() function

}  // End of class
