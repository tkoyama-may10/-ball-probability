## Synopsis

This program computes the ball probability and the Fisher-Bingham integral by the holonomic gradient method (HGM).
For details, see `note-en.tex`.

## Code Example

In /bin/bash/, the following will compute the value of the ball probability 
for dim=3, r=20.0, tau=(3.0, 2.0, 1.0) and lambda=(1.0, 0.5, 0.25).

`./a.out 1 3 20.0 3.0 2.0 1.0 1.0 0.5 0.25`

Here, the first argument switches the behavior of the program.
To compute the value of the Fisher-Bingham integral at the same parameters,
change the first argument of the above command into '10'.

## Motivation

On this project, we release programs and raw data used in our paper [1].

## Installation

Download the zip file `-ball-probability-master.zip`.
Then run the following commands:
```
unzip ball-probability-master.zip
cd ball-probability-master/
latex note-en.tex && latex note-en.tex
gcc program.c -lm -lgsl -lblas -O0 -Wall 
```
For compiling `program.c`, you will need the GNU scientific library and BLAS.


## API Reference

1. T.Koyama, A.Takemura, 
``Holonomic gradient method for distribution function of 
a weighted sum of noncentral chi-square random variables'',
http://arxiv.org/abs/1503.00378.
2. GNU scientific library, http://www.gnu.org/software/gsl/

## Tests

```
$ ./a.out 1 3 2.5 3.0 2.0 1.0 1.0 0.5 0.25
Probability=   0.319980
$ ./a.out 10 3 2.5 3.0 2.0 1.0 1.0 0.5 0.25
r=1.000000, C=4.18879e-06
f:		   96045.8     105034     181305     866601     842933     733832 
f:		   90855.4    99357.8     171507     819770     797381     694176 
r=2.500000
f:		   82334.6     191619     171096     301573 1.56444e+06     721234 
g:		   3.44882    8.02651    7.16684    12.6322    65.5312     30.211 
lap:		   30.7793    20.9859    8.24445    30.7793    182.196    35.1861 
ratio:		   0.11205   0.382472   0.869292   0.410414   0.359675   0.858605 
Fisher-Bingham integral:	27.951000
  Laprace approximation:	30.779291
(FB)/(Laprace):			0.908111
```


## Contributors

Tamio Koyama

## License

GPL
