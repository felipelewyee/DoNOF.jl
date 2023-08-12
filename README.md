# DoNOF.jl

This is a Julia version of the [DoNOF](https://github.com/DoNOF/DoNOFsw/) (Donostia Natural Orbital Functionals) software written by Prof. Mario Piris, with the intention of take advantage of the wonderful capabilities of Julia language.

# Requirements

You should have Julia installed. 

Example:
~~~
wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.2-linux-x86_64.tar.gz
tar zxvf julia-1.9.2-linux-x86_64.tar.gz
sudo mv julia-1.9.2 /opt
~~~
Append to .bashrc
~~~
export PATH="$PATH:/opt/julia-1.9.2/bin"
~~~
or simple do: apt install julia

# Installation (for internal development)

For development, you only need to clone DoNOF.jl from github and change to the project directory
~~~
git clone https://github.com/felipelewyee/DoNOF.jl.git
cd DoNOF.jl
~~~
[Suggestion] Install Revise.jl (only the first time)
~~~
julia
]
activate .
add Revise
~~~
To develop DoNOF.jl
~~~
julia
]
activate .
# Return to green prompt (press backspace key)
using Revise # Assuming that it is installed
using DoNOF
~~~
Example of input
~~~
using DoNOF

mol = """
0 1
  O  0.0000   0.000   0.121
  H  0.0000   0.751  -0.485
  H  0.0000  -0.751  -0.485
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=true)

p.ipnof = 8

p.RI = true
p.gpu = false

p.method = "ID"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true,do_mbpt=false,do_translate_to_donofsw=true)
~~~

# Installation (for normal execution)

Open a julia prompt and do:
~~~
add https://github.com/felipelewyee/DoNOF.jl # Only the first time
using DoNOF
~~~
then you can do
~~~
input("example.jl")
~~~
or simply copy and paste the content of the file in the prompt.
You can also do
~~~
julia example.jl > example.out
~~~
