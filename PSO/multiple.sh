#!/bin/bash

$(javac -cp contest.jar player51.java)
$(jar cmf MainClass.txt submission.jar player51.class)

for i in {400..500};
do
  echo $(java -jar testrun.jar -submission=player51 -evaluation=KatsuuraEvaluation -seed=$i) >> results_katsuura.txt
done
