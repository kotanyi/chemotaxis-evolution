#!/bin/bash
#$1 - the output file to be parsed
#opens the file specified by $1, reads the protein ids and creates an R script to plot the file (otherwise i would have to manually change the protein ids in a generic R script).

echo 'args <- commandArgs(TRUE)
outputs <- read.csv(file = args[1], head = TRUE, sep = ",")
jpeg_filename <- paste(args[1], ".jpg", sep = "", collapse = NULL)
jpeg(filename = jpeg_filename, width = 1500, height = 1500, quality = 75, res = 200)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 6))
plot(outputs$Time, outputs$Protein4, type = "l", col = "green")
par(mar = c(5, 4, 4, 2) + 0.1)' > plot_protein_activities.R
counter=1

for i in $(head -n 1 $1 | grep -o 'Protein[0-9]\+' | grep -vw 'Protein4'); do
	echo 'lines(outputs$Time, outputs$'$i', col = "red", lty = '$counter')' >> plot_protein_activities.R
	counter=$((counter + 1))
done

echo 'time_range <- range(outputs$Time)' >> plot_protein_activities.R
echo -n 'legend(time_range[2] + 1, 0.5, c(' >> plot_protein_activities.R

for i in $(head -n 1 $1 | grep -o 'Protein[0-9]\+' | grep -vw 'Protein4'); do
	protein_id=$(echo $i | grep -o '[0-9]\+')
	echo -n '"'$protein_id'", ' >> plot_protein_activities.R
done

echo -n '"4", "input"), col = c(' >> plot_protein_activities.R
counter=1

while [ $counter -le $(head -n 1 $1 | grep -o 'Protein[0-9]\+' | grep -vw 'Protein4' | wc -l) ]; do
	echo -n '"red", ' >> plot_protein_activities.R
	counter=$((counter + 1))
done

echo -n '"green", "orange"), lty = c(' >> plot_protein_activities.R
counter=1

while [ $counter -le $(head -n 1 $1 | grep -o 'Protein[0-9]\+' | grep -vw 'Protein4' | wc -l) ]; do
	echo -n $counter', ' >> plot_protein_activities.R
	counter=$((counter + 1))
done

echo '1, 1))
dev.off()' >> plot_protein_activities.R
