# Coding_Examples
Coding Examples from college and other projects

RMS2D.m
Matlab program to calculate 2D spatial RMS using a variety of weighting functions.
Initially used to calculate fracture intensity from curcature attributes extracted from seismic, but can be used on any spatial data.
Imports tabular data using the XYZimporter.m program below.

FRACAZ.m
Matlab program which takes spatial data and plots the orientations of any identifiable features within the data.
Originally used to produce rose diagrams of fracture orientations from curvature data extracted from seismic.
Could be used to identify any linear feature orientation within spatial data, e.g.
Dominant road orientation within a city, grain structure of wood using photographs, etc.
Imports tabular data using the XYZimporter.m program below.

XYZimporter.m
Imports tabular data using Matlab interfaces. Data for import should be in XYZ (also known as XYAttribute) format.
Output is a matrix of the data with NaN fill for blank spaces, and a structure field with pertinent grid information.


XYZexporter.m
Exports matrix data from matlab to tabular XYZ (also known as XYAttribute) format.
other required info comes from the structured field portion of the XYZimporter.


GeophoneNoiseFilter.m
Matlab program from my Master's thesis. Removes air noise from geophone recorded seismic signals using a microphone.


CRAZYALG.m
Work in progress...
designed to solve systems of linear equations where the variables are not set scalars, rather they are probablility distributions.
M=a*X+b*Y+c*Z
N=d*X+e*Y+f*Z
O=g*X+h*Y+i*Z

Where M, N, and O are known outputs.
Where a through i are known coefficients
Where X, Y, and Z are the unknown probability distributions, each with a mean and standard deviation.

An example problem to be solved with this program is this:
take 1 trip to 3 different grocery stores and blindly buy (i.e. without looking at individual prices) an amount of apples, bananas, and oranges at each store.
the total price paid at each store are the outputs (M, N, O).
The number of each fruit purchased at each store are the coefficients (a trhrough i).
X, Y, and Z are the prices for each of the fruits.
However, the prices will change between stores so X, Y, and Z are not set numbers and cannot be solved using traditional methods.
Rather X, Y, and Z are distributiuons of prices. each with their own mean and standard deviation
CRAZYALG.m solves for each mean and standard deviation using a modified genetic algorithm and monte carlo simulation.
