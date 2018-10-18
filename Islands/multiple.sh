#!/bin/bash

$(javac -cp contest.jar player51.java)

for i in {0..1};
do
  echo $(java -jar testrun.jar -submission=player51 -evaluation=SchaffersEvaluation -seed=$i >~/Desktop/EvoCompetent/Islands/results_schaffers/results_seed$i.txt)
done
