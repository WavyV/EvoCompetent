#!/bin/bash

$(javac -cp contest.jar player51.java)
$(jar cmf MainClass.txt submission.jar player51.class)

for i in {10..99};
do
  echo $(java -jar testrun.jar -submission=player51 -evaluation=KatsuuraEvaluation -seed=$i >~/Desktop/ec_assignment/Unimodal/results_katsuura/results_seed$i.txt)
done
