
mkdir -p resources/homer
cd resources/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install hg38 mm10