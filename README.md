## Synopsis

This program computes the ball probability and the Fisher-Bingham integral 
by the holonomic gradient method (HGM).
For details, see `note-en.tex`.

## Code Example

In /bin/bash/, the following will compute the value of the ball probability 
for dim=3, r=20.0, sigma=diag(3.0, 2.0, 1.0) and mu=(1.0, 0.5, 0.25).

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
$ ./a.out 3 3 2.5 3.0 2.0 1.0 1.0 0.5 0.25
f:		   8.04362    7.48803    6.68603    73.6548    61.1349    28.1842 
lap:		   19.4788    25.8904    8.09076    179.795    230.425    34.6386 
ratio:		  0.412943    0.28922   0.826378   0.409661   0.265314   0.813665 
Fisher-Bingham integral:	26.0758
  Laprace approximation:	28.7672
(FB)/(Laprace):			0.906444
```


## Contributors

Tamio Koyama

## License

GPL
