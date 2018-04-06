- This is a repo that I made in response to
  [this question on Bioinformatics Stack Exchange - What is the fastest way to get the reverse complement of a DNA sequence in python?](https://bioinformatics.stackexchange.com/questions/3583)
- To run this test on your own machine, navigate to this directory and execute `python reverse_complement_tests.py`
- Feel free to contribute if you have other ideas. Exceedingly
  convoluted or bad implementations for fun are also welcome.
  
- To date, here is a table of the speeds when run on my Apple machine.

```
      the runtime of reverse complement implementations.
    10000 strings and 250 repetitions
    ╔══════════════════════════════════════════════════════╗
    ║                            %inc   s total  str per s ║
    ╠══════════════════════════════════════════════════════╣
    ║ name                                                 ║
    ║ user172818 seqpy.c        93.7%  0.002344  4266961.4 ║
    ║ alexreynolds Cython (v2)  93.4%  0.002468  4051583.1 ║
    ║ alexreynolds Cython (v1)  90.4%  0.003596  2780512.1 ║
    ║ devonryan string          86.1%  0.005204  1921515.6 ║
    ║ jackaidley bytes          84.7%  0.005716  1749622.2 ║
    ║ jackaidley bytesstring    83.0%  0.006352  1574240.6 ║
    ║ biopython just rc         70.5%  0.012275   814638.6 ║
    ║ global dict                5.4%  0.035330   283046.7 ║
    ║ revcomp_translateSO       45.9%  0.020202   494999.4 ║
    ║ string_replace            37.5%  0.023345   428364.9 ║
    ║ revcom from SO            28.0%  0.026904   371694.5 ║
    ║ naive (baseline)           1.5%  0.036804   271711.5 ║
    ║ lambda from SO           -39.9%  0.052246   191401.3 ║
    ║ biopython seq then rc    -32.0%  0.049293   202869.7 ║
    ╚══════════════════════════════════════════════════════╝
```

# Installation

- Install [pipenv](https://docs.pipenv.org/)
- run `pipenv install`
