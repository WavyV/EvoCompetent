#!/bin/bash

$(javac -cp contest.jar player51.java)

for i in {0..99};
do
  echo $(java -jar testrun.jar -submission=player51 -evaluation=BentCigarFunction -seed=$i >~/Desktop/ec_assignment/PSO/results_bentcigar/results_seed$i.txt)
done
