function fchk(filename,p,bset,jobtype,E_t,elag,n,C)

    ne = 0
    for i in 1:bset.natoms
      ne -= bset.atoms[i].Z
    end
    charge = ne + p.nalpha + p.nbeta

    max_l = 0
    max_nprimitives = 0
    nprimitives = 0
    for i in 1:bset.nshells
       max_l = max(max_l,bset.basis[i].l)
       max_nprimitives = max(max_nprimitives,size(bset.basis[i].exp)[1])
       nprimitives += size(bset.basis[i].exp)[1]
    end

    shells2atom = zeros(bset.nshells)
    i = 0
    for iatom in 1:bset.natoms
        nshells_atom = bset.shells_per_atom[iatom]
        shells2atom[i+1:i+nshells_atom] .= iatom
        i += nshells_atom
    end

    f = open(filename*".fchk","w")

    @printf(f,"%s\n",filename)
    if(p.ista==0)
        @printf(f,"%s PNOF%d %s\n",jobtype,p.ipnof,bset.name)
    else
        @printf(f,"%s PNOF%ds %s\n",jobtype,p.ipnof,bset.name)
    end
    @printf(f,"Number of atoms                            I           %6d\n",p.natoms)
    @printf(f,"Charge                                     I           %6d\n",charge)
    @printf(f,"Multiplicity                               I           %6d\n",p.mul)
    @printf(f,"Number of electrons                        I           %6d\n",p.ne)
    @printf(f,"Number of alpha electrons                  I           %6d\n",p.nalpha)
    @printf(f,"Number of beta electrons                   I           %6d\n",p.nbeta)
    @printf(f,"Number of basis functions                  I           %6d\n",p.nbf)
    @printf(f,"Number of independant functions            I           %6d\n",p.nbf5)
    @printf(f,"Number of contracted shells                I           %6d\n",bset.nshells)
    @printf(f,"Highest angular momentum                   I           %6d\n",max_l)
    @printf(f,"Largest degree of contraction              I           %6d\n",max_nprimitives)
    @printf(f,"Number of primitive shells                 I           %6d\n",nprimitives)
#    print("Virial Ratio                               R                {}".format(p.natoms),file=f)
#    print("SCF Energy                                 R                {}".format(p.natoms),file=f)
    @printf(f,"Atomic numbers                             I   N=      %6d\n",bset.natoms)
    for i in 1:bset.natoms
        Z = bset.atoms[i].Z
        @printf(f," %11d",Z)
        if((i)%6==0 || i==p.natoms)
            @printf(f,"\n")
	end
    end
    @printf(f,"Nuclear Charges                            R   N=      %6d\n",bset.natoms)
    for i in 1:bset.natoms
        Z = bset.atoms[i].Z
        @printf(f," %.8e",Z)
        if((i)%6==0 || i==p.natoms)
            @printf(f,"\n")
	end
    end
    @printf(f,"Current cartesian coordinates              R   N=      %6d\n",bset.natoms*3)
    idata = 0
    for i in 1:bset.natoms
        for ixyz in 1:3
            idata += 1
            @printf(f," % .8e",bset.atoms[i].xyz[ixyz]*1.88973)
            if(idata%5==0 || idata==p.natoms*3)
                @printf(f,"\n")
            end
        end
    end
    @printf(f,"Shell types                                I   N=      %6d\n",bset.nshells)
    for ishell in 1:bset.nshells
        if(p.spherical && bset.basis[ishell].l > 1)
	    @printf(f," %11d",-1*bset.basis[ishell].l)
        else
            @printf(f," %11d",bset.basis[ishell].l)
	end
	if((ishell)%6==0 || ishell==bset.nshells)
            @printf(f,"\n")
        end
    end
    @printf(f,"Number of primitives per shell             I   N=      %6d\n",bset.nshells)
    for ishell in 1:bset.nshells
        @printf(f," %11d",size(bset.basis[ishell].exp)[1])
        if((ishell)%6==0 || ishell==bset.nshells)
            @printf(f,"\n")
	end
    end
    @printf(f,"Shell to atom map                          I   N=      %6d\n",bset.nshells)
    for ishell in 1:bset.nshells
        @printf(f," %11d",shells2atom[ishell])
	if((ishell)%6==0 || ishell==bset.nshells)
            @printf(f,"\n")
        end
    end

    @printf(f,"Coordinates of each shell                  R   N=      %6d\n",bset.nshells*3)
    idata = 0
    for ishell in 1:bset.nshells
        for ixyz in 1:3
            idata += 1
            @printf(f," % 8e",bset.basis[ishell].atom.xyz[ixyz])
            if(idata%5==0 || idata==bset.nshells*3)
                @printf(f,"\n")
	    end
        end
    end

    @printf(f,"Total Energy                               R     %.15e\n",E_t)
    @printf(f,"Primitive exponents                        R   N=      %6d\n",nprimitives)
    idata = 0
    for ishell in 1:bset.nshells
	for iprim in 1:size(bset.basis[ishell].exp)[1]
            idata += 1
	    @printf(f," % .8e",bset.basis[ishell].exp[iprim])
            if(idata%5==0 || idata==nprimitives)
                @printf(f,"\n")
            end
        end
    end
    @printf(f,"Contraction coefficients                   R   N=      %6d\n",nprimitives)
    idata = 0
    for ishell in 1:bset.nshells
	for iprim in 1:size(bset.basis[ishell].exp)[1]
            idata += 1
            @printf(f," % .8e",bset.basis[ishell].coef[iprim])
            if(idata%5==0 || idata==nprimitives)
                @printf(f,"\n")
            end
        end
    end

    e_val = diag(elag)
    @printf(f,"Alpha Orbital Energies                     R   N=      %6d\n",p.nbf)
    for i in 1:p.nbf
        @printf(f," % .8e",e_val[i])
        if((i)%5==0 || i==p.nbf)
            @printf(f,"\n")
        end
    end
    @printf(f,"Alpha MO coefficients                      R   N=      %6d\n",p.nbf*p.nbf)    
    idata = 0
    for j in 1:p.nbf
        for i in 1:p.nbf
            idata += 1
            @printf(f," % .8e",C[i,j])
            if(idata%5==0 || idata==p.nbf*p.nbf)
                @printf(f,"\n")
	    end
        end
    end
    @printf(f,"Total SCF Density                          R   N=       %6d\n",p.nbf*(p.nbf+1)/2)
    Cnbf5 = C[1:end,1:p.nbf5]
    @tullio DM[m,nn] := 2*n[i]*Cnbf5[m,i]*Cnbf5[nn,i]
    idata = 0
    for mu in 1:p.nbf
        for nu in 1:mu
            idata += 1
	    @printf(f," % .8e",DM[mu,nu])
            if(idata%5==0 || idata==p.nbf*(p.nbf+1)/2)
                @print(f,"\n")
	    end
        end
    end
    @printf(f,"Natural orbital occupancies                R   N=       %6d\n",p.nbf5)
    for i in 1:p.nbf5
        @printf(f," % .8e",n[i])
        if(i%5==0 || i==p.nbf5)
            @printf(f,"\n")
        end
    end

    close(f)

end
