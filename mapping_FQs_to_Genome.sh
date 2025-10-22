##PULL PIPELINE
nextflow pull Edward-ward-762/postPlasmidsaurusQC -r main

##RUN PIPELINE
nextflow run Edward-ward-762/postPlasmidsaurusQC \
-r main \
-profile docker \
-resume \
--inputFile ./inputFile.csv