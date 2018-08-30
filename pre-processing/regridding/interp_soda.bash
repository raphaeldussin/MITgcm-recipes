# first run one year to create weights
./interp_SODA3.3.1_ASTE.py 1981
# then run the rest in parallel processes
seq -f %4g 1982 2015 | parallel -j 24 ./interp_SODA3.3.1_ASTE.py {}
