#!/bin/bash

$(javac -cp contest.jar player51.java)
$(jar cmf MainClass.txt submission.jar player51.class)

for i in {1..100};
do
  echo $(java -jar testrun.jar -submission=player51 -evaluation=SchaffersEvaluation -seed=$i) >> results_schaffers.txt
done
