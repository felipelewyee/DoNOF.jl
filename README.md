# <img src="https://github.com/felipelewyee/DoNOF.jl/blob/master/DoNOFjl.png" height=100>

This is a Julia version of the [DoNOF](https://github.com/DoNOF/DoNOFsw/) (Donostia Natural Orbital Functionals) software written by Prof. Mario Piris, with the intention of take advantage of the wonderful capabilities of Julia language.

## 🌟 Installation

Open a julia prompt, enter to the Pkg REPL by pressing ], and simply add DoNOF:
~~~julia
add DoNOF
~~~

## 🎯 Example

- Put the following directly on the Julia prompt to run a calculation:
~~~julia
using DoNOF

mol = """
0 1
 O  0.0000   0.000   0.121
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

bset,p = DoNOF.molecule(mol,"cc-pvtz")
p.ipnof = 8

DoNOF.energy(bset,p)
~~~

- To save the input in a file and run it directly on the command line:
~~~bash
julia example.jl > example.out
~~~

## 😎 Tips

- You should have [Julia installed](https://julialang.org/downloads)

- Remember to configure the number of threads for parallelism, for example:
~~~
export JULIA_NUM_THREADS=12
~~~

- To use GPUs, just add CUDA, cuTENSOR and KernelAbstractions:
~~~julia
using CUDA, cuTENSOR, KernelAbstractions, DoNOF

mol = """
0 1
 O  0.0000   0.000   0.121
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

bset,p = DoNOF.molecule(mol,"cc-pvtz")
p.ipnof = 8

DoNOF.energy(bset,p)
~~~

