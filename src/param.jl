#export Param

mutable struct Param
    natoms::Int64
    nbf::Int64
    nbfaux::Int64
    nalpha::Int64
    nbeta::Int64
    ne::Int64
    mul::Int64
    ndoc::Int64
    nsoc::Int64
    ndns::Int64
    nvir::Int64
    closed::Bool
    nac::Int64
    nbf5::Int64
    no0::Int64
    title::String
    maxit::Int64
    thresheid::Float64
    no1::Int64
    maxitid::Int64
    maxloop::Int64
    ipnof::Int64
    ista::Int64
    threshl::Float64
    threshe::Float64
    threshec::Float64
    threshen::Float64
    scaling::Bool
    nzeros::Int64
    nzerosm::Int64
    nzerosr::Int64
    itziter::Int64
    diis::Bool
    thdiis::Float64
    ndiis::Int64
    perdiis
    ncwo::Int64
    noptorb::Int64
    nv::Int64
    gradient::String
    optimizer::String
    gpu::Bool
    RI::Bool
    HighSpin::Bool
    MSpin::Int64
    lamb::Float64
    
end

function Param(bas_name)

    Fermi.Options.set("basis", bas_name)

    bset = Fermi.GaussianBasis.BasisSet()

    natoms = bset.natoms
    nbf = bset.nbas
    nalpha = bset.molecule.Nα
    nbeta = bset.molecule.Nβ
    mul = bset.molecule.multiplicity
   
    nbfaux = 0	
    ne = nalpha + nbeta
    no1 = 0

    for i in 1:natoms
	Z = bset.molecule.atoms[i].Z
        if 1<=Z && Z<=  2
            no1 += 0           # H-He
        elseif 3<=Z && Z<= 10
            no1 +=  1          # Li-Ne
        elseif 11<=Z && Z<= 18
            no1 +=  5          # Na-Ar
        elseif 19<=Z && Z<= 36
            no1 +=  9          # K-Kr
        elseif 37<=Z && Z<= 49
            no1 += 18          # Rb-In
        elseif 50<=Z && Z<= 54
            no1 += 23          # Sn-Xe
        elseif 55<=Z && Z<= 71
            no1 += 27          # Cs-Lu
        elseif 72<=Z && Z<= 81
            no1 += 30          # Hf-Tl
        elseif 82<=Z && Z<= 86
            no1 += 39          # Pb-Rn
        elseif 87<=Z && Z<=109
            no1 += 43          # Fr-Mt
        end
    end

    ndoc = nbeta - no1
    nsoc = nalpha - nbeta
    ndns = ndoc + nsoc
    nvir = nbf - nalpha

    ncwo = 1
    if ne == 2
        ncwo = -1
    end
    if ndns!=0
        if ndoc>0
	    if ncwo!=1
	        if ncwo==1 && ncwo>nvir/ndoc
		    ncwo = trunc(Int, nvir/ndoc)
                end
	    end
	else
	    ncwo = 0
        end
    end

    closed = (nbeta == (ne+mul-1)/2 && nalpha == (ne-mul+1)/2)

    noptorb = nbf
    nac = ndoc * (1+ncwo)
    nbf5 = no1 + nac + nsoc
    no0 = nbf - nbf5

    title = "pydonof"
    maxit = 10000  # Número máximo de iteraciones de Occ-SCF
    thresheid = 10^-6#8 # Convergencia de la energía total
    maxitid = 300  # Número máximo de iteraciones externas en HF
    maxloop = 30  # Iteraciones internas en optimización orbital
    ipnof = 7     # PNOFi a calcular
    ista = 0     # PNOFi a calcular
    threshl = 10^-5#4   # Convergencia de los multiplicadores de Lagrange
    threshe = 10^-6#6   # Convergencia de la energía
    threshec = 10^-10 # Convergencia  de la energía en optimización orbital
    threshen = 10^-10 # Convergencia  de la energía en optimización de ocupaciones
    scaling = true     # Scaling for f
    nzeros = 0
    nzerosm = 5
    nzerosr = 0
    itziter = 10        # Iteraciones para scaling constante
    diis = true         # DIIS en optimización orbital
    thdiis = 10^-3     # Para iniciar DIIS
    ndiis = 5           # Número de ciclos para interpolar matriz de Fock generalizada en DIIS
    perdiis = true      # Aplica DIIS cada NDIIS (true) o después de NDIIS (false)
    ncwo = ncwo         # Número de orbitales débilmente ocupados acoplados a cada orbital fueremtente ocupado
    noptorb = noptorb   # Número de orbitales a optimizar Nbf5 <= Noptorb <= Nbf
    nv = ncwo*ndoc
    gradient = "analytical"
    optimizer = "BFGS"
    gpu = false
    RI = false

    HighSpin = false
    MSpin = 0

    lamb = 0.0

    return Param(natoms,
    nbf,
    nbfaux,
    nalpha,
    nbeta,
    ne,
    mul,
    ndoc,
    nsoc,
    ndns,
    nvir,
    closed,
    nac,
    nbf5,
    no0,
    title,
    maxit,
    thresheid,
    no1,
    maxitid,
    maxloop,
    ipnof,
    ista, 
    threshl, 
    threshe, 
    threshec, 
    threshen, 
    scaling, 
    nzeros, 
    nzerosm, 
    nzerosr, 
    itziter, 
    diis, 
    thdiis, 
    ndiis, 
    perdiis, 
    ncwo, 
    noptorb, 
    nv, 
    gradient, 
    optimizer, 
    gpu,
    RI,
    HighSpin,
    MSpin,
    lamb)
end

function set_ncwo(p,ncwo)
    if p.ne==2
        ncwo= -1
    end
    if p.ndns!=0
        if p.ndoc>0
            if ncwo!=1
                if ncwo==-1 || ncwo > p.nvir/p.ndoc
                    ncwo = trunc(Int,p.nvir/p.ndoc)
		end
	    end
        else
            ncwo = 0
	end
    end

    p.ncwo = ncwo

    p.nac = p.ndoc * (1 + ncwo)
    p.nbf5 = p.no1 + p.nac + p.nsoc   #JFHLY warning: nbf must be >nbf5
    p.no0 = p.nbf - p.nbf5
    p.nv = p.ncwo*p.ndoc

end
