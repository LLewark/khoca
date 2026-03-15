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

## 1. Installation


### 1.1 Download of binaries

Binaries for Linux are available for download from 
http://lewark.de/lukas/khoca.html
They should run on any Linux installation that has python3.6. Binaries for
Windows or Mac are not available at the moment.

### 1.2 Run in a Docker container

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

### 1.3 Installation from [PyPI](https://pypi.org/project/khoca/) via `pip`

This installation method provides an interactive access to Khoca inside a Python session or program.

```bash
pip install khoca
```

or for a specific version (1.4 in the example):

```bash
pip install khoca==1.4
```

This also works with [SageMath](https://www.sagemath.org/):

```bash
sage -pip install khoca
```

### 1.4 Source code

The source code, including instructions on how to compile it, is available at
the GitHub repository khoca:
https://github.com/LLewark/khoca

## 2. Usage

### 2.1 Usage at the `bash` prompt

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

### 2.2 Usage inside a Python session or program

To use Khoca inside a Python session, function or method of a Python class you need to install it via `pip` (see install hint 1.3 above). In Python you can import an interactive Khoca calculator.

It's output consists of two lists of quadruples the first for reduced and the second for unreduced homology. In such a quadruple the first item stands for the `t`-degree, the second for the `q`-degree, the third for the torsion and the last item stands for the coefficient of the corresponding summand of the homology.

If the option `print_messages` is given, than the output contains a second result which is a list of all print messages of the command-line version.

Examples:

```Python console
    >>> from khoca import InteractiveCalculator
    >>> KH = InteractiveCalculator()
    >>> KH
    Khovanov homology calculator for Frobenius algebra: Z[X] / (1*X^2).
    >>> KH('braidaBaB')
    [[[-2, 4, 0, 1], [-1, 2, 0, 1], [0, 0, 0, 1], [1, -2, 0, 1], [2, -4, 0, 1],
      [-2, 4, 0, 0], [-1, 4, 0, 0], [-1, 2, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0],
      [1, 0, 0, 0], [1, -2, 0, 0], [2, -2, 0, 0]], [[-2, 3, 0, 1], [-2, 5, 0, 1],
      [-1, 1, 0, 1], [-1, 3, 0, 1], [0, -1, 0, 1], [0, 1, 0, 1], [1, -3, 0, 1],
      [1, -1, 0, 1], [2, -5, 0, 1], [2, -3, 0, 1], [-1, 3, 2, 1], [-2, 3, 0, -1],
      [-1, 3, 0, -1], [-2, 5, 0, 0], [-1, 5, 0, 0], [-1, 1, 0, 0], [0, 1, 0, 0],
      [-1, 3, 0, 0], [0, 3, 0, 0], [0, -1, 0, 0], [1, -1, 0, 0], [0, 1, 0, 0],
      [1, 1, 0, 0], [2, -3, 2, 1], [1, -3, 0, -1], [2, -3, 0, -1], [1, -1, 0, 0],
      [2, -1, 0, 0]]]
    >>> res, mess = KH('braidaaa', print_messages=True); mess # doctest: +NORMALIZE_WHITESPACE
    ['Result:', 'Reduced Homology:', 'Non-equivariant homology:',
     't^-3q^8 + t^-2q^6 + t^0q^2', 'Unreduced Homology:',
     'Non-equivariant homology:', 't^-3q^9 + t^-2q^5 + t^0q^1 + t^0q^3 + t^-2q^7[2]']
    >>> KH((1, 1, 1)) == KH('braidaaa')
    True
    >>> KH([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]) == KH('braidaBaB')
    True
    >>> KH((1, -2, 1, -2)) == KH('braidaBaB')
    True
    >>>
```

In these examples default values for the base ring, the frobenius algebra, the root and equivariance are used. For special values of these properties you need to define specific instances of the interactive calculator. In an `IPython` session you may obtain more information by online-help using `?`:

```Python console
In [1]: from khoca import InteractiveCalculator

In [2]: InteractiveCalculator?


Init signature:
InteractiveCalculator(
    coefficient_ring=0,
    frobenius_algebra=(0, 0),
    root=0,
    equivariant=None,
)
Docstring:     
Class to allow the usage of `Khoca` interactively in a Python session.

EXAMPLES::

    >>> from khoca import InteractiveCalculator
... (as above)
    >>> KH((1, -2, 1, -2)) == KH('braidaBaB')
    True
    >>>
Init docstring:
Constructor.

INPUT:

    - ``coefficient_ring`` -- coefficient ring ``F``of the homology. It
      can be given by an integer (or string convertible to an integer):
      -  ``0`` (default) the ring of integers
      -  ``1``  the rational field
      -  a prime for the corresponding finite field
    - ``frobenius_algebra`` -- the Frobenius algebra ``F[x]/p`` given by
      the coefficients of the normed polynomial ``p`` (default is
      ``p = X^2``) as a tuple of length ``deg(p)`` where the constant
      term is the first entry.
    - ``root`` -- a root of the polynomial ``p`` (default is zero).
    - ``equivariant`` optional integer ``n > 1`` giving the degree of
      ``sl(n)`` for equivariant homology. If this is given then the input
      to the keyword arguments  ``frobenius_algebra`` and ``root`` will
      be ignored. If it is not given then non equivariant homology will
      be calculated.

EXAMPLES::

    >>> from khoca import InteractiveCalculator
    >>> InteractiveCalculator(1, (0, 1), 0)
    Khovanov homology calculator for Frobenius algebra: Q[X] / (1*X^2 + 1*X).
    >>> from khoca import InteractiveCalculator
    >>> KH = InteractiveCalculator(equivariant=3); KH
    Khovanov homology calculator for Frobenius algebra: Z[a, b, c][X] / (1*X^3 + c*X^2 + b*X + a).
```

To obtain the documentation of the instance-call apply `?` to an instance of the calculator:

```Python console
In [3]: KH = InteractiveCalculator()

In [4]: KH?
Signature:      KH(link, command=None, print_messages=False, verbose=False, progress=False)
Type:           InteractiveCalculator
...  (as above)
    Khovanov homology calculator for Frobenius algebra: Z[a, b, c][X] / (1*X^3 + c*X^2 + b*X + a).
Call docstring:
Instance call to apply the calculator to a link with respect to a
certain command.

INPUT:

    - ``link`` -- the link to which the calculation should be done.
      It can be given as a Tuple which will be interpreted as a braid
      in Tietze form which closure is the link in question. Further
      a list of lists is accepted which will be interpreted as a list
      of crossings in pd-notation. Alternatively you can declare the
      link by a string. The following formats are accepted:

      - ``BraidX`` for a braid formatted as in knotscape (``a`` = first
        Artin generator, ``A`` = its inverse, ``b`` = second Artin
        generator, etc.). This works only for ``sl(2)`` homology,
        otherwise output is nonsensical
      - ``PdX`` for input in PD notation (as e.g. given on !KnotInfo).
        Again, this works only for ``sl(2)`` homology, otherwise output
        is nonsensical
      - ``GaussX`` for bipartite knot given as a matched diagram,
        following the convention explained in the section below. This
        works for ``sl(N)`` homology for all ``N``
      - ``command`` -- a command given as string as explained in the
        command line version of the functionality. This can be:
      - ``MirrorX`` takes the dual of the result at spot ``X``
      - ``SumXY`` computes the homology of the connected sum of the
        results saved at spots ``X`` and ``Y`` (numbers separated by
        a non-digit character)
      - ``CalcX` outputs the result saved at spot X. If you forget
        this command, the program will have no output

      If this keyword argument is not given it will be interpreted
      as ``Calc0`` by default.

      - ``print_messages`` boolean (default is ``False``). If set to
        ``True`` the print output of the command line version is
        returned as a list on the second output position. By default
        all print messages to ``stdout`` of the command line version
        are suppressed
      - ``verbose`` boolean (default is ``False``). If it is set to
        ``True`` all print messages to ``stdout`` together with special
        verbose messages of the command line version are printed
      - ``progress`` boolean (default is ``False``). If it is set to
        ``True`` progress bars will be printed

OUTPUT:

    Two lists of quadruples the first for reduced and the second for
    unreduced homology. In such a quadruple the first item stands for
...
EXAMPLES::

    >>> from khoca import InteractiveCalculator
    >>> KH = InteractiveCalculator(1, (1, 0), 0); KH
    Khovanov homology calculator for Frobenius algebra: Q[X] / (1*X^2 + 1).
    >>> KH('braidaaa', 'dual0 sum0+1 braidaBaB sum2+3 calc4', print_messages=True)
    ([[[-5, 10, 0, 1], [-4, 8, 0, 1], [-4, 8, 0, 1], [-3, 6, 0, 1],
       [-3, 6, 0, 1], [-3, 6, 0, 1], [-2, 4, 0, 1], [-2, 4, 0, 1],
       [-2, 4, 0, 1], [-2, 4, 0, 1], [-2, 4, 0, 1], [-2, 4, 0, 1],
       [-1, 2, 0, 1], [-1, 2, 0, 1], [-1, 2, 0, 1], [-1, 2, 0, 1],
       [-1, 2, 0, 1], [-1, 2, 0, 1], [-1, 2, 0, 1], [0, 0, 0, 1],
       [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 1],
       [0, 0, 0, 1], [0, 0, 0, 1], [1, -2, 0, 1], [1, -2, 0, 1],
       [1, -2, 0, 1], [1, -2, 0, 1], [1, -2, 0, 1], [1, -2, 0, 1],
       [1, -2, 0, 1], [2, -4, 0, 1], [2, -4, 0, 1], [2, -4, 0, 1],
       [2, -4, 0, 1], [2, -4, 0, 1], [2, -4, 0, 1], [3, -6, 0, 1],
       [3, -6, 0, 1], [3, -6, 0, 1], [4, -8, 0, 1], [4, -8, 0, 1],
       [5, -10, 0, 1]], [[0, -1, 0, 1], [0, 1, 0, 1]]],
     ['Result:', 'Reduced Homology:', 'Non-equivariant homology:',
      'Page 1:', 't^-5q^10 + 2t^-4q^8 + 3t^-3q^6 + 6t^-2q^4 + 7t^-1q^2 + 7t^0q^0
               + 7t^1q^-2 + 6t^2q^-4 + 3t^3q^-6 + 2t^4q^-8 + t^5q^-10',
      'The spectral sequence collapses on the first page.\n',
      'Unreduced Homology:', 'Non-equivariant homology:',
      'Page 1:', 't^-5q^11 + t^-4q^7 + t^-4q^9 + t^-3q^5 + 2t^-3q^7 + 2t^-2q^3
               + 4t^-2q^5 + 4t^-1q^1 + 3t^-1q^3 + 4t^0q^-1 + 4t^0q^1 + 3t^1q^-3
               + 4t^1q^-1 + 4t^2q^-5 + 2t^2q^-3 + 2t^3q^-7 + t^3q^-5 + t^4q^-9
               + t^4q^-7 + t^5q^-11',
      'Page 3 = infinity:', 't^0q^-1 + t^0q^1'])
```

## 3. Encoding of matched diagrams

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

## 4. References

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
