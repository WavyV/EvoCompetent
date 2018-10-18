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
		} else {
			// Do sth else
		}
	}
	
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
		} // End for loop
		
		individual = population[indexInd];
		
		return individual;
	} // End fittestEnd function
	
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
		} // End for loop
		
		individual = population[indexInd];
		
		return individual;
	} // End LeastFit function
	
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
	} //End indexMax function
	
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


	public void run()
	{
		int evals = 0;
		int N = 180;  // Number of individuals, must be divisible by 12 
		int D = 10;  // Dimensions
		int NrImmigrants = 2;
		int ImmigrationGen = 18; // After this number of generations immigration takes place
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

		
		// Method-specific parameters for Unimodal
		double sigma_init = 2.5;  // Sigma of Gaussian used in initialization
		double tau = 0.3;  // Sigma of Gaussian used in mutation
		double epsilon1 = 0.00001;  // Minimum value for sigma
		
		// Method-specific parameters for Differential Evolution
		double F = 0.2;  // Scaling factor
		double Cr = 0.5;  // Crossover probability
		
		// Method-specific parameters for PSO
		boolean multikulti = true; // True means we use the base method, False means we use the consensus method
		double w = 0.8;  // inertia
		double phi1 = 0.1;  // learning rate for personal influence
		double phi2 = 0.1;  // learning rate for social influence
		double epsilon3 = 0.1;
		
		
		

		
		
		//Initializing population1 (Unimodal)
		
		for (int i=0; i<(N/3); i++) {
			for (int j=0; j<D; j++) {
				population1[i][j] = rnd_.nextGaussian()*sigma_init;
			}
			for (int k=D; k<(3*D)+1; k++) {
				population1[i][k]=0;
			}
			population1[i][(3*D)+1] = sigma_init;
			population1[i][(3*D)+2] = 0;
			population1[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(population1[i], 0, D));
			evals++;
		}
		
		//Initializing population2 (DE)
		
		for (int i=0; i<(N/3); i++) {
			for (int j=0; j<D; j++) {
				population2[i][j] = rnd_.nextDouble()*5;
				if (rnd_.nextDouble()<0.5) {
					population2[i][j] = -1 * population2[i][j];
				}
			}
			//for (int k=D; k<(3*D)+3; k++) {
			//	population2[i][k]=0;
			//}
			population2[i][(3*D)+3] = (double) evaluation_.evaluate(Arrays.copyOfRange(population1[i], 0, D));
			evals++;				
		}
		
		//Initializing population3 (PSO)
		
		for (int i=0; i<(N/3); i++) {
			for (int j=0; j<D; j++) {
				population3[i][j] = rnd_.nextDouble()*5;
				if (rnd_.nextDouble()<0.5) {
					population3[i][j] = -1 * population3[i][j]; 
				}
				population3[i][(2*D)+j] = population3[i][j];
			}
			for (int k=D; k<(2*D); k++) {
				population3[i][k] = rnd_.nextGaussian()*3;
			}
			for (int h=(2*D); h<(3*D)+3; h++) {
				population3[i][h] = 0;
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
		
		
		
		
		// Start evolution
		
		while (evals + N < evaluations_limit_) { //+180 because otherwise we run out of evaluations during the while loop, which gives an error
			
			
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
				children[i][D] = children[i][D] * Math.exp(tau*rnd_.nextGaussian());
				if(children[i][D] < epsilon1){
					children[i][D] = epsilon1;
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
			
			if(evals >= evaluations_limit_-1){
				System.out.println(Arrays.toString(GenerationChampion));
			}
			
			
		     
		} //end of while loop
	} //end of run function
} //end of class
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		/*

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
		
		*/

		/*
		while(evals < evaluations_limit_) {
			
			if (evals % 25 =! 1) {

				for(int i=0; i<N; i++){
	
					double [] perturbed_velocity = new double[D];
					for(int j=0; j<D; j++){
						perturbed_velocity[j] = w * individuals[i][j+D] + phi1 * rnd_.nextDouble() * (individuals[i][j+2*D] - individuals[i][j]) + phi2 * rnd_.nextDouble() * (champion[j] - individuals[i][j]);
					}
	
					for(int j=0; j<D; j++){
						individuals[i][j] += perturbed_velocity[j];
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
			}
			else {
				//make sure new immigrant at PSO doesn't change
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

		}  // End of while loop

		
		
		
		
		int NrImmigrants = 2;
		int ImmigrationEvals = 25; // Every 25 evaluations immigration takes place
		double [][] immigrants = new double [NrImmigrants * NrIslands][D+1];
		double [] ranges = new double [10];
		//double [][] population1 = new double [N][D];
		//double [][] population2 = new double [N][D];
		//double [][] population3 = new double [N][D];
		Random r = new Random();
		Random s = new Random();
		Random t = new Random();
		Random u = new Random();
		Random v = new Random();
		Random w = new Random();
		
		*/
		
		/*
		
		if (evals % ImmigrationEvals == 0) { 
			if (multikulti == true) {
				immigrants[0]=fittestInd(population1)[0];
				int index1 = (int) LeastFit(population1)[1][0];
				immigrants[2]=fittestInd(population2)[0];
				int index2 = (int) LeastFit(population2)[1][0];
				immigrants[4]=fittestInd(population3)[0];
				int index3 = (int) LeastFit(population3)[1][0];
				
				int firstindex = r.nextInt(34);
				int secondindex = s.nextInt(34);
				int thirdindex = t.nextInt(35);

				immigrants[1]=population1[firstindex];
				immigrants[3]=population2[secondindex];
				immigrants[5]=population3[thirdindex];
				
				
				for (int h=immigrants[5].length-1 ; h>D-1; h--) {
					immigrants[4][h]=0;
					immigrants[5][h]=0;
				}
				
				int fourthindex = u.nextInt(34);
				int fifthindex = v.nextInt(34);
				int sixthindex = w.nextInt(35);
				
				population2[index2]=immigrants[0]; 
				population2[fourthindex]=immigrants[1];
				population3[index3]=immigrants[2];
				population3[fifthindex]=immigrants[3];
				population1[index1]=immigrants[4];	
				population1[fourthindex]=immigrants[5];
					
				}
			else {
				for (int i=0; i<D; i++) {
					double [] counting = new double [10];
					for (int n=0; n<10; n++) {
						counting[n]=0;
					}
					for (int k=0; k<D; k++) {
						for (int l=0; l<10; l++) {
							if (population1[k][i] < (-4+l) && population1[k][i] >= (-5+l)) {
								counting[l]++;
							}
						}
					}
					ranges[i]=indexMax(counting);
					
					
				}
					
				}
						
			}
			
	*/	

	//}  // End of run() function
	
	

//}  // End of class
