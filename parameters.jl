module parameters
export Param

mutable struct Param
    natoms
    nbf
    nbfaux
    nalpha
    nbeta
    ne
    mul
    ndoc
    nsoc
    ndns
    nvir
    nac
    nbf5
    no0
    title
    maxit
    thresheid
    no1
    maxitid
    maxloop
    ipnof
    ista
    threshl
    threshe
    threshec
    threshen
    scaling
    nzeros
    nzerosm
    nzerosr
    itziter
    diis
    thdiis
    ndiis
    perdiis
    ncwo
    noptorb
    nv
    gradient
    optimizer
    gpu
    RI
    HighSpin
    MSpin
    
end

function Param(natoms,nbf,nalpha,nbeta,mul,Z_list)
   
    nbfaux = 0	
    ne = nalpha + nbeta
    no1 = 0
    for i in 1:natoms
	Z = Z_list[i]
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
		    ncwo = int(nvir/ndoc)
                end
	    end
	else
	    ncwo = 0
        end
    end

    noptorb = nbf
    #closed = nbeta ==
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
    nzeros = 1
    nzerosm = 5
    nzerosr = 2
    itziter = 10        # Iteraciones para scaling constante
    diis = true         # DIIS en optimización orbital
    thdiis = 10^-3     # Para iniciar DIIS
    ndiis = 5           # Número de ciclos para interpolar matriz de Fock generalizada en DIIS
    perdiis = true      # Aplica DIIS cada NDIIS (true) o después de NDIIS (false)
    ncwo = ncwo         # Número de orbitales débilmente ocupados acoplados a cada orbital fueremtente ocupado
    noptorb = noptorb   # Número de orbitales a optimizar Nbf5 <= Noptorb <= Nbf
    nv = ncwo*ndoc
    gradient = "analytical"
    optimizer = "CG"
    gpu = false
    RI = false

    HighSpin = false
    MSpin = 0
    

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
    MSpin)
end




end



