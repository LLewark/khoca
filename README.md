# README for Khoca, a knot homology calculator

The same text as in this README file is available (with illustrations) at
http://lewark.de/lukas/khoca.html

Khoca is computer program to calculate sl(N)-homology of knots. The program has
been written for joint projects with Andrew Lobb such as [3, 4]. The paper [3]
also contains a description of the algorithm used by khoca. The main innovation
is to use Krasner's calculation of the sl(N)-homology of the basic two-crossing
tangle [2] for calculations of the homology of bipartite knots.

Khoca calculates the following:

* Khovanov sl(2)-homology of arbitrary links, given as a braid or in PD code.
* Khovanov-Rozansky sl(N)-homology with N > 2 of bipartite knots, given by a
  certain encoding of a matched diagram of the knot (see [3] and section
  "Encoding of matched diagrams" below).
* Homology over the integers, the rationals or a prime field.
* Either equivariant homology, or homology with an arbitrary fixed potential.
* All pages of the spectral sequence of filtered homology over a field.
* Reduced and unreduced homology.
* Homology of sums and mirror images of knots.

You are encouraged to contact me with any kind of questions or comments
regarding khoca. If you are using khoca for a project or publication, please
cite this web page, or the paper [3].

## Installation

Binaries for Linux are available for download from 
http://lewark.de/lukas/khoca.html
They should run on any Linux installation that has python3.6. Binaries for
Windows or Mac are not available at the moment.

The source code, including instructions on how to compile it, is available at
the GitHub repository khoca:
https://github.com/LLewark/khoca

To run Khoca in Docker type:

```bash
docker run -it soehms/khoca:latest
```

Its download size is 162 MB and it will need 516 MB of disk space on your
device. To create a new (resp. locally own) Docker image cd to the khoca
directory type

```bash
docker build -f Dockerfile --tag khoca:<your_tag> .
```

If your machine has an older CPU it can happen that you get *Illegal Instruction*
errors. In that case you better should use the image `soehms/khoca:old_cpu`.


## Usage

To use the program, run khoca.py (a python3 script) from the command line.
khoca.py takes three arguments:

1. The coefficient ring; `0` for integers, `1` for rationals, a prime `p` for the
corresponding finite field.

2. A sequence of N integers `a_0, ..., a_{N-1}` separated by a non-digit
character, defining the Frobenius algebra `F[X]/(X^N + a_{N-1}X^{N-1} + ... +
a_0)`.  Alternatively, `e` followed by a number `N` for equivariant computation
over `sl(N)`. For example, `-1.0.0` gives the Frobenius algebra `F[X]/(X3 - 1)`.

3. A root of the polynomial given in 2. for the calculation of reduced homology
(for the dependence of reduced homology on a root, see [3]).  For example, to
get the standard graded reduced homology, use 0 as root. If you are not
interested in reduced homology, it does not matter what root you chose (and
khoca does not check that the number is actually a root).

The option `-p` will show progress bars, `-v` will give more verbose
non-mathematical information, and `-h` will print a short help text. Each
argument after the first three arguments, can be one of the following.

1. `BraidX` calculates homology of a link given as closure of the braid `X`,
formatted as in knotscape (`a` = first Artin generator, `A` = its inverse, `b` =
second Artin generator, etc.). This works only for `sl(2)` homology, otherwise
output is nonsensical.

2. `PdX` calculates homology of a link given in PD notation (as e.g. given on
KnotInfo). Again, this works only for `sl(2)` homology, otherwise output is
nonsensical.

3. `GaussX` calculates homology of a bipartite knot given as a matched diagram,
following the convention explained in the section below. This works for `sl(N)`
homology for all `N`.

4. `MirrorX` takes the dual of the result at spot `X`.

5. `SumXY` computes the homology of the connected sum of the results saved at
spots `X` and `Y` (numbers separated by a non-digit character).

6. `CalcX` outputs the result saved at spot `X`. If you forget this command, the
program will have no output.

The program keeps a stack of computed homologies, enumerated 0,1,2... . Each of
the commands 1 - 5 puts a new homology on that stack, whereas the command 6.
prints the homology at a certain spot. This is mainly useful to compute
homology of sums of knots.

Here are some examples:

```bash
./khoca.py 0 0.0 0 braidaBaB calc0
```

calculates the classical `sl(2)` Khovanov homology (both reduced and unreduced)
of the closure of the braid `aBaB` (knotscape notation), i.e. the figure-eight
knot.


```bash
./khoca.py 0 e2 0 pd[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]] calc0
```

calculates integral equivariant `sl(2)` homology of the figure-eight knot.


```bash
./khoca.py 7 0.-1 0 braidabcdefabcdefabcdefabcdefabcdefabcdefabcdefabcdef
calc0 -p
```

calculates Khovanov homology of the `(7,8)`-torus knot over `F7` with perturbed
potential, displaying progress bars. This calculation takes roughly two
minutes, and shows that the spectral sequence does not collapse on the second
page, refuting the knight-move conjecture over finite fields (cf. [1]).


```bash
./converters/montesinos.py [1/5,1/3,-1/2]
```

outputs `[12,4,16,10,15,9,14,13]`, the code for a matched diagram of the
`(5,3,-2)`-pretzel knot, aka the `(3,5)`-torus knot, aka `10_{124}`. So

```bash
./khoca.py 1 1.0.0.0.0 0 gauss[12,4,16,10,15,9,14,13] calc0 calculates
```

rational `sl(5)` homology and the corresponding Rasmussen invariant of the `(3,5)`
torus knot.


```bash
./khoca.py 1 1.0 0 braidaaa dual0 sum0+1 braidaBaB sum2+3 calc4
```

calculates `sl(2)` homology of the sum of the trefoil, its mirror image and a figure-8-knot.


# Encoding of matched diagrams

This section describes how to encode a matched knot diagram, i.e. a diagram
that consists of `n` copies of the basic 2-crossing tangle. Resolving each basic
tangle into two intervals and a chord results in a single circle with n
non-intersecting (red) chords, which may be on either side of the circle.
Enumerate the `2n` chord endpoints by walking around the circle. If a chord
connects the points `i` and `j`, let `f(i) = j`. Write down the list `f(1), f(2), ...,
f(2n)` omitting `f(i)` if `f(i) < i`. Moreover, make the list entries signed, and
let the sign reflect the sign of the two crossings of the corresponding
2-crossing tangle. This list of `n` non-zero integers uniquely determines the
matched diagram. As an example, the standard diagram of `6_1` is matched and
encoded as `[-4,6,5]`.

"Half" of Montesinos knots are bipartite [3]. You may use the python3 script
`./converters/montesinos.py` to obtain the encoding of a matched diagram of
Montesinos knots.

# References

[1] Bar-Natan: Fast Khovanov Homology Computations, Journal of Knot Theory and
its Ramifications 16 (2007), no.3, pp. 243-255, arXiv:math/0606318, MR2320156.

[2] Daniel Krasner: A computation in Khovanov-Rozansky Homology, Fundamenta
Mathematicae 203 (2009), pp. 75-95, arXiv:0801.4018, MR2491784.

[3] Lukas Lewark and Andrew Lobb: New Quantum Obstructions to Sliceness,
Proceedings of the London Mathematical Society 112 (2016), no. 1, pp. 81-114,
arXiv:1501.07138, MR3458146.

[4] Lukas Lewark and Andrew Lobb: Upsilon-like concordance invariants from
sl(n) knot cohomology, arXiv:1707.00891.


Lukas Lewark, 2018
