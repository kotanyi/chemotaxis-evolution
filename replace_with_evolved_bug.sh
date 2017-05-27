#!/bin/bash
#extract a bug from the last generation of an effector_motor_rates file and paste the information into a parameters file
#$1 - specifies a effector_motor_rates1_?_pre file to be used

cp parameters1 parameters1_new
new_value=$(tail -n 1 $1 | awk -F ',' '{ print $2 }')
sed -i 's/effectorMotorAssociationRate = .\+/effectorMotorAssociationRate = '$new_value'/' parameters1_new
new_value=$(tail -n 1 $1 | awk -F ',' '{ print $3 }')
sed -i 's/effectorMotorDissociationRate = .\+/effectorMotorDissociationRate = '$new_value'/' parameters1_new
new_value=$(tail -n 1 $1 | awk -F ',' '{ print $4 }')
sed -i 's/realisticBugConstant = .\+/realisticBugConstant = '$new_value'/' parameters1_new
new_value=$(tail -n 1 $1 | awk -F ',' '{ print $5 }')
sed -i 's/realisticBugLinear = .\+/realisticBugLinear = '$new_value'/' parameters1_new
new_value=$(tail -n 1 $1 | awk -F ',' '{ print $6 }')
sed -i 's/realisticBugQuadratic = .\+/realisticBugQuadratic = '$new_value'/' parameters1_new
new_value=$(tail -n 1 $1 | awk -F ',' '{ print $7 }')
sed -i 's/realisticBugMemoryLength = .\+/realisticBugMemoryLength = '$new_value'/' parameters1_new
