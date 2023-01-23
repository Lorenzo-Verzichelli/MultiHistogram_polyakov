# MultiHistogram_polyakov

This code reweights data from simulation to extrapolate or interpolate the results at other temperatures

In particular I'm intrested in the Polyakov Loop and its susceptivity in a Yang Mills lactice simulation

The code is based on a program I recived from Lorenzo Maio which solved a similar task

MultiHistogram.cpp and .h file contain the code I recived. I have done very little change
The theoretical principles can be found in chap. 8 of Newman Barkema

Jackknife.cpp and .h file contain some function usefull to estimate the error on resampled observables
I wrote these imitating Claudio Bonati's python code for the same task
(avaiable on Git Hub: https://github.com/claudio-bonati/statanalysis/blob/master/multihistogram.py)

weighted_mean.cpp and .h contain the implementation of a class used in the jackknife analysis

to compile I have a bash script outside the repository folder with the following instruction:

cd MultiHistogram_polyakov
g++ -O3 -Wall -Wextra -Werror -pedantic -Wconversion -o ../MultiHisto Jackknife.cpp MultiHistogram.cpp main.cpp
