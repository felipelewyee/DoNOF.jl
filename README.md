# DoNOF.jl

This is a Julia version of the [DoNOF](https://github.com/DoNOF/DoNOFsw/) (Donostia Natural Orbital Functionals) software written by Prof. Mario Piris, with the intention of take advantage of the wonderful capabilities of Julia language.

# <img src="https://github.com/felipelewyee/DoNOF.jl/blob/master/DoNOFjl.png" height=150>

# Requirements

You should have [Julia installed](https://julialang.org/downloads)

# Installation (for normal execution)

Open a julia prompt, enter to the Pkg REPL by pressing ], and simply add DoNOF:
~~~
add DoNOF
~~~

# Example

You can copy-paste the following directly on the Julia prompt to run a calculation:
~~~
using DoNOF

mol = """
0 1
 O  0.0000   0.000   0.121
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

bset,p = DoNOF.molecule(mol,"cc-pvtz")

p.ipnof = 8

p.RI = true
p.gpu = false

p.orb_method = "Rotations"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p)
~~~

You can also save the input in a file and run it directly on the command line:
~~~
julia example.jl > example.out
~~~

Remember to have configured the number of threads that Jullia can use for parallelism, for example:
~~~
export JULIA_NUM_THREADS=12
~~~
