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
    double [][] individuals = new double[N][3*D+2];  // 2D matrix used to keep track of individuals

		// STRUCTURE OF INDIVIDUALS:
		// First 10 entries: Indices 0 u/i 9 -> Point in space which is evaluated
		// Second 10 entries: Indices 10 u/i 19 -> Velocity vector
		// Third 10 entries: Indices 20 u/i 29 -> The coordinates of the personal best point in space
		// Index 30: the fitness of the personal best
		// Index 31: the fitness of the current point in space

    // Method-specific parameters
    double w = 0.4;  // inertia
    double phi1 = 0.9;  // learning rate for personal influence
    double phi2 = 1.2;  // learning rate for social influence
		double epsilon = 0.001;
    double [] champion = new double[D+1];


    // Initialize matrix randomly and initial fitness evaluation
    for(int i=0; i<N; i++){
      for(int j=0; j<D; j++){
        individuals[i][j] = rnd_.nextDouble()*5;
        individuals[i][j+2*D] = individuals[i][j];
      }
      for(int j=D; j<2*D; j++){
        individuals[i][j] = rnd_.nextGaussian()*3;
      }
      individuals[i][3*D] = (double) evaluation_.evaluate(Arrays.copyOfRange(individuals[i], 0, D));  
      individuals[i][3*D+1] = individuals[i][3*D];
			evals++;
    }

    // Find first champion
    int index_best = 0;
    for(int i=0; i<N; i++){
      if(individuals[i][3*D] > individuals[index_best][3*D]){
        index_best = i;
      }
    }
    for(int i=0; i<D; i++){
      champion[i] = individuals[index_best][i];
    }
    champion[D] = individuals[index_best][3*D];

		// if(champion[D] < initial[D]){
		// 	champion = initial;
		// }


    while(evals < evaluations_limit_) {

      for(int i=0; i<N; i++){

        double [] perturbed_velocity = new double[D];
        for(int j=0; j<D; j++){
          perturbed_velocity[j] = w * individuals[i][j+D] + phi1 * rnd_.nextDouble() * (individuals[i][j+2*D] - individuals[i][j]) + phi2 * rnd_.nextDouble() * (champion[j] - individuals[i][j]);
					if(perturbed_velocity[j] <= epsilon){
						if(perturbed_velocity[j] < 0){
							perturbed_velocity[j] = -epsilon;
						} else {
							perturbed_velocity[j] = epsilon;
						}
					}

				}

        for(int j=0; j<D; j++){
          individuals[i][j] += perturbed_velocity[j];
					if(individuals[i][j] < -5){
						individuals[i][j] = -5;
					}
					if(individuals[i][j] > 5){
						individuals[i][j] = 5;
					}
        }

        for(int j=D; j<2*D; j++){
          individuals[i][j] = perturbed_velocity[j-D];
        }
        individuals[i][3*D+1] = (double) evaluation_.evaluate(Arrays.copyOfRange(individuals[i], 0, D));
				evals++;
				if(individuals[i][3*D+1] > individuals[i][3*D]){
          individuals[i][3*D] = individuals[i][3*D+1];
          for(int j=2*D; j<3*D; j++){
            individuals[i][j] = individuals[i][j-2*D];
          }
        }

      }

      // Find champion
      for(int i=0; i<N; i++){
        if(individuals[i][3*D] > champion[D]){
          for(int j=0; j<D; j++){
            champion[j] = individuals[i][j];
          }
          champion[D] = individuals[i][3*D];
        }
      }

			if(evals >= evaluations_limit_-1){
				System.out.println(Arrays.toString(champion));
			}

    }  // End of while loop

	}  // End of run() function

}  // End of class
