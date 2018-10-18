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
		// Set seed of algortihm"s random process
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
	}


	// This function will find and return the index of some value in an array (assuming it only occurs once)
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


	// This function ranks the population from worst (index 0) to best (index N-1) based on their fitness
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


	// What does this do?
	public static double[] fittestInd(double[][] population)
	{

		int popSize = population.length;
		int lengthArray = population[0].length;
		double [] individual = new double[lengthArray];
		double fitness = population[0][lengthArray-1];
		int indexInd = 0;

		for (int j=1; j<popSize; j++) {
			if (population[j][lengthArray-1]>fitness) {
				fitness = population[j][lengthArray-1];
				indexInd = j;
			}
		}

		individual = population[indexInd];

		return individual;
	}



	// What does this do?
	public static double[] LeastFit(double[][] population)
	{
		int size = population.length;
		int lengthArray = population[0].length;
		double [] individual = new double[lengthArray];
		double fitness = population[0][lengthArray-1];
		int indexInd = 0;

		for (int j=0; j < size; j++) {
			if (population[j][lengthArray-1]<fitness) {
				fitness = population[j][lengthArray-1];
				indexInd = j;
			}
		}

		individual = population[indexInd];

		return individual;
	}


	// What does this do?
	public static int LeastFitIndex(double[][] population) {
		int D = 10;
		double lowestFit = population[0][(3*D)+3];
		int len = population.length;
		int index = 0;
		for (int i=1; i<len; i++) {
			if (population[i][(3*D)+3]<lowestFit) {
				index = i;
			}
		}
		return index;
	}


	// What does this do?
	public static int indexMax(double [] IntList) {
		double maxi = IntList[0];
		int maxIndex = 0;
		int len = IntList.length;
		for (int i=1; i<len; i++) {
			if (IntList[i]>maxi) {
				maxi=IntList[i];
				maxIndex=i;
			}
		}
		return maxIndex;
	}


	// What does this do?
	public static double sigmaAverage(double [][] population) {
		int D = 10;
		double total = 0;
		int len = population.length;
		for (int i=0; i<len; i++) {
			total += population[i][(3*D)+1];
		}
		double average = total/len;

		return average;
	}


	// Method to print a matrix to the console
	public static void printArray(double matrix[][])
	{
		for (double[] row : matrix)
			System.out.println(Arrays.toString(row));
	}




	// Main method in which the algorithm is actually run
	public void run()
	{
		int evals = 0;  // Current number of fitness evaluations done
		int N = 180;  // Number of individuals, must be divisible by 12
		int D = 10;  // Dimensions
		int NrImmigrants = 2;
		int ImmigrationGen = 20; // After this number of generations immigration takes place
		int GenCount = 0; //Count number of generations to make sure immigration takes place after every ... generations
		int imm = 0;
		double [][] population1 = new double[N/3][(3*D)+4];  // 2D matrix used to keep track of individuals in population 1
		double [][] population2 = new double[N/3][(3*D)+4];  // 2D matrix used to keep track of individuals
		double [][] population3 = new double[N/3][(3*D)+4];  // 2D matrix used to keep track of individuals

		/*
		Index 0-(D-1) are the values of the 10 dimensions
		Index D-(2D-1) is the perturbated velocity
		Index 2D-(3D-1) are the values of the dimensions for the personal best fitness
		Index 3D gives the fitness of the personal best position
		Index 3D+1 gives the sigma value of the individual
		Index 3D+2 tells whether the individual is a recent immigrant (1) or not (0)
		Index 3D+3 is the fitness of the individual at the current position
		 */


		// Method-specific parameters for Conventional Algorithm
		double tau = 0.3;  // Sigma of Gaussian used in mutation
		double epsilon1 = 0.0001;  // Minimum value for sigma

		// Method-specific parameters for Differential Evolution
		double F = 0.2;  // Scaling factor
		double Cr = 0.5;  // Crossover probability
		double epsilon2 = 0.0001;  // Minimum value for perturbation

		// Method-specific parameters for PSO
		boolean multikulti = true; // True means we use the base method, False means we use the consensus method
		double w = 0.5;  // inertia
		double phi1 = 1.5;  // learning rate for personal influence
		double phi2 = 1.2;  // learning rate for social influence
		double epsilon3 = 0.0001;  // Minimum value for velocity




		//Initializing population1 (CA)
		for (int i=0; i<(N/3); i++) {
			for (int j=0; j<D; j++) {
				population1[i][j] = -5 + rnd_.nextDouble()*10;  // Locations
			}
			population1[i][(3*D)+1] = rnd_.nextDouble()*2.5;  // Sigmas
			population1[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(population1[i], 0, D));
			evals++;
		}

		//Initializing population2 (DE)
		for (int i=0; i<(N/3); i++) {
			for (int j=0; j<D; j++) {
				population2[i][j] = -5 + rnd_.nextDouble()*10;  // Locations
			}
			population2[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(population1[i], 0, D));
			evals++;
		}

		//Initializing population3 (PSO)
		for (int i=0; i<(N/3); i++) {
			for (int j=0; j<D; j++) {
				population3[i][j] = -5 + rnd_.nextDouble()*10;  // Locations
				population3[i][(2*D)+j] = population3[i][j];  // Personal best
			}
			for (int k=D; k<(2*D); k++) {
				population3[i][k] = -2.5 + rnd_.nextDouble()*5;  // Velocities
			}
			population3[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(population1[i], 0, D));
			population3[i][3*D] = population3[i][(3*D)+3];
			evals++;
		}


		//Lists to keep track of the over all champion and champion per generation
		double [] GenerationChampion1 = new double [(3*D)+4];
		double [] GenerationChampion2 = new double [(3*D)+4];
		double [] GenerationChampion3 = new double [(3*D)+4];
		double [] GenerationChampion = new double [(3*D)+4];
		double fitnesschampion;
		double [] Alltime_champion3 = fittestInd(population3);


		System.out.println("===============NewGeneration");
		System.out.println(GenCount);
		System.out.println("===============ConventionalAlgorithm");
		printArray(population1);
		System.out.println("===============DifferentialEvolution");
		printArray(population2);
		System.out.println("===============ParticleSwarmOptimisation");
		printArray(population3);


		// Start evolution. This is the main evolutionary cycle. The magic happens here.
		while (evals + N < evaluations_limit_) { //evals+180 because otherwise we run out of evaluations during the while loop, which gives an error


			if (GenCount % ImmigrationGen == 0 && GenCount != 0) {

				//immigration
				//evals stays the same

				double [][] immigrants = new double [NrImmigrants * 3][(3*D)+4];
				int [] indices = new int [NrImmigrants*3];
				double aveSigma = sigmaAverage(population1);
				imm = 1;

				for (int i=0; i<(N/3); i++) { //Setting immigration checker to zero for everyone
					population1[i][(3*D)+2]=0;
					population2[i][(3*D)+2]=0;
					population3[i][(3*D)+2]=0;
				}



				//First choosing the fittest individual of every population to be copied


				immigrants[0]=fittestInd(population1);
				immigrants[NrImmigrants]=fittestInd(population2);
				immigrants[2*NrImmigrants]=fittestInd(population3);

				//And calculating the indices of the least fit individuals they should replace

				indices[0] = LeastFitIndex(population2);
				indices[NrImmigrants] = LeastFitIndex(population3);
				indices[2*NrImmigrants] = LeastFitIndex(population1);


				//Then choosing other random individuals to be copied and adjusted or replaced

				for (int j=0; j<3; j++) { // Nr of islands

					for (int i=1; i<NrImmigrants; i++) {

						if (j==0) { //Immigrants from island 1 (Uni) to island 2 (DE)
							immigrants[i] = population1[(int) rnd_.nextDouble()*(N/3)];
							indices[i] = (int) rnd_.nextDouble()*(N/3);
						}
						if (j==1) { //Immigrants from island 2 (DE) to island 3 (PSO)
							immigrants[NrImmigrants+i] = population2[(int) rnd_.nextDouble()*(N/3)];
							for (int k=0; k<D; k++) {
								immigrants[NrImmigrants+i][D+k] = 0; // Potential old velocity perturbation is set to zero
								immigrants[NrImmigrants+i][(2*D)+k] = immigrants[2+i][k]; // Personal best position is current position
							}
							immigrants[NrImmigrants+i][3*D] = immigrants[2+i][(3*D)+3]; // Personal best fitness is current fitness
							indices[NrImmigrants+i] = (int) rnd_.nextDouble()*(N/3);
						}
						if (j==2) { //Immigrants from island 3 (PSO) to island 1 (Uni)
							immigrants[(2*NrImmigrants)+i] = population3[(int) rnd_.nextDouble()*(N/3)];
							immigrants[(2*NrImmigrants)+i][(3*D)+1] = aveSigma; // Setting sigma to average of population
							indices[(2*NrImmigrants)+i] = (int) rnd_.nextDouble()*(N/3);
						}
					}
				}

				for (int i=0; i<(NrImmigrants*3); i++) {

					immigrants[i][(3*D)+2]=1; //Setting immigration checker to one for every immigrant

				}

				// And finally actually immigrating

				for (int k=0; k<NrImmigrants; k++) {
					population2[indices[k]] = immigrants[k];
					population3[indices[NrImmigrants+k]] = immigrants[NrImmigrants+k];
					population1[indices[(2*NrImmigrants)+k]] = immigrants[(2*NrImmigrants)+k];
				}



			} //End of if statement


			// The evolution of population 1 (Uni)

			population1 = rank_population(population1);

			int parent_number = 1;
			int [] dummy_array = new int [N/3];
			double [] parent1 = new double [(3*D)+4];
			double [] parent2 = new double [(3*D)+4];
			double [][] children = new double [(N/3)/2][(3*D)+4];

			for (int i=0; i<((N/3)/2); i++) { //We pick (N/3)/2 parents in total

				int random_index = 0;
				boolean picked = false;

				while (picked == false) {

					random_index = rnd_.nextInt((N/3)/2) + ((N/3)/2);
					if (dummy_array[random_index] == 0) {

						dummy_array[random_index] = 1;
						picked = true;

					}
				}

				if (parent_number == 1){
			        parent1 = population1[random_index];
			        parent_number++;
			    } else {
			        parent2 = population1[random_index];
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
				children[i][3*D+1] = children[i][3*D+1] * Math.exp(tau*rnd_.nextGaussian());
				if(children[i][3*D+1] < epsilon1){
					children[i][3*D+1] = epsilon1;
				}
			}
			// Apply mutation, simple fixed sigma Gaussian mutation with 50% probability on each gene
		     for(int i=0; i<children.length; i++){
		       for(int j=0; j<D; j++){
		         if(rnd_.nextDouble() < 1){
		           children[i][j] += rnd_.nextGaussian()*children[i][3*D+1];
		         }
		       }
		     }

		      // Check fitness of children
		     for(int i=0; i<children.length; i++){
		       children[i][3*D+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(children[i], 0, D));
					evals++;
		     }

		     // Select survivors
				// First put all the previous generation and the children together in one matrix
				// Then rank the matrix based on fitness
				// The select the best N to form the new generation
		     double [][] total_population = new double[population1.length+children.length][D+2];
		     for(int i=0; i<population1.length; i++){
		    	 total_population[i] = population1[i];
		     }
		     for(int i=0; i<children.length; i++){
		    	 total_population[i+population1.length] = children[i];
		     }
		     total_population = rank_population(total_population);

				// Selecting the new generation at the end of the array to have it ranked properly
		     int index = population1.length-1;
		     for(int i=total_population.length-1; i >= children.length; i--){
		    	 population1[index] = total_population[i];
		    	 index--;
		     } // End of evolution population1 (Uni)


		     // The evolution of population2 (DE)

		     double [][] mutant_population = new double [N/3][(3*D)+4];
		     double [][] trial_population = new double [N/3][(3*D)+4];



		     for(int i=0; i<(N/3); i++){
		    	 double [] base_vector = population2[rnd_.nextInt(N/3)];
		    	 double [] parent_1 = population2[rnd_.nextInt(N/3)];
		    	 double [] parent_2 = population2[rnd_.nextInt(N/3)];

		    	 double [] perturbation_vector = new double[(3*D)+4];
		    	 for(int j=0; j<D; j++){
					perturbation_vector[j] = F * (parent_1[j] - parent_2[j]);
					mutant_population[i][j] = base_vector[j] + perturbation_vector[j];
		    	 }

		    	 for(int j=0; j<D; j++){
		    		 if(rnd_.nextDouble() < Cr){
		    			 trial_population[i][j] = mutant_population[i][j];
		    		 } else {
		    			 trial_population[i][j] = population2[i][j];
		    		 }
		    	 }

		    	 int randIndex = rnd_.nextInt(D);
		    	 trial_population[i][randIndex] = parent_1[randIndex];
		     }



		     for(int i=0; i<(N/3); i++){
				trial_population[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(trial_population[i], 0, D));
				evals++;
				if(trial_population[i][(3*D)+3] > population2[i][(3*D)+3]){
					population2[i] = trial_population[i];
				}
		     } // End of evolution population2 (DE)


		     //The evolution of population3 (PSO)



		     for (int i=0; i<(N/3); i++) {
		    	 double [] perturbed_velocity = new double[D];
		    	 for (int j=0; j<D; j++) {

		    		 if (imm == 1) { //Keeps track if immigration just happened (imm=1) or not (imm=0)

		    			 if (population3[i][(3*D)+2] == 0) {

		    				 perturbed_velocity[j] = w * population3[i][j+D] + phi1 * rnd_.nextDouble() * (population3[i][j+2*D] - population3[i][j]) + phi2 * rnd_.nextDouble() * (Alltime_champion3[j] - population3[i][j]);
		    				 if (Math.abs(perturbed_velocity[j])<=epsilon3) {
		    					 if (perturbed_velocity[j]<0) {
		    						 perturbed_velocity[j] = -1 * epsilon3;
		    					 } else {
		    						 perturbed_velocity[j] = epsilon3;
		    					 }
		    				 }

		    			 }

		    		 } else {

		    			 perturbed_velocity[j] = w * population3[i][j+D] + phi1 * rnd_.nextDouble() * (population3[i][j+2*D] - population3[i][j]) + phi2 * rnd_.nextDouble() * (Alltime_champion3[j] - population3[i][j]);
	    				 if (Math.abs(perturbed_velocity[j])<=epsilon3) {
	    					 if (perturbed_velocity[j]<0) {
	    						 perturbed_velocity[j] = -1 * epsilon3;
	    					 } else {
	    						 perturbed_velocity[j] = epsilon3;
	    					 }
	    				 }

		    		 }

		    	 }

		    	 for (int j=0; j<D; j++) {
		    		 population3[i][j] += perturbed_velocity[j];
		    		 if (population3[i][j] < -5) {
		    			 population3[i][j] = -5;
		    		 }
		    		 if (population3[i][j] > 5) {
		    			 population3[i][j] = 5;
		    		 }
		    		 population3[i][D+j] = perturbed_velocity[j];
		    	 }
		    	 population3[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(population3[i], 0, D));
		    	 evals++;
		    	 if (population3[i][(3*D)+3] > population3[i][3*D]) {
		    		 population3[i][3*D] = population3[i][(3*D)+3];
		    		 for (int k=0; k<D; k++) {
		    			 population3[i][(2*D)+k]=population3[i][k];
		    		 }
		    	 }

		     }

		     if (fittestInd(population3)[(3*D)+3]>Alltime_champion3[(3*D)+3]) {
		    	 Alltime_champion3 = fittestInd(population3);
		     }


			// Updating general variables and champions

			imm = 0;
			GenCount++;

			GenerationChampion1 = fittestInd(population1);
			GenerationChampion2 = fittestInd(population2);
			GenerationChampion3 = fittestInd(population3);
			fitnesschampion = 0;
			if (GenerationChampion1[(3*D)+3] > fitnesschampion) {
				fitnesschampion = GenerationChampion1[(3*D)+3];
				GenerationChampion = GenerationChampion1;
			}
			if (GenerationChampion2[(3*D)+3] > fitnesschampion) {
				fitnesschampion = GenerationChampion2[(3*D)+3];
				GenerationChampion = GenerationChampion2;
			}
			if (GenerationChampion3[(3*D)+3] > fitnesschampion) {
				fitnesschampion = GenerationChampion3[(3*D)+3];
				GenerationChampion = GenerationChampion3;
			}

			System.out.println("===============NewGeneration");
			System.out.println(GenCount);
			System.out.println("===============ConventionalAlgorithm");
			printArray(population1);
			System.out.println("===============DifferentialEvolution");
			printArray(population2);
			System.out.println("===============ParticleSwarmOptimisation");
			printArray(population3);



		} //end of while loop
	} //end of run function
} //end of class
