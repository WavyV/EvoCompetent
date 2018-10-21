#!/bin/bash

$(javac -cp contest.jar player51.java)

for i in {0..99};
do
  echo $(java -jar testrun.jar -submission=player51 -evaluation=SchaffersEvaluation -seed=$i >~/documents/GitHub/EvoCompetent/Islands/Results_Schaffers/results_seed$i.txt)
done
