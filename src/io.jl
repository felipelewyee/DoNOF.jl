export molecule

function read_xyz_from_file(filename)

    f = open(filename,"r")

    filecontent = read(f,String)
    filecontent = split(filecontent,"\n")

    natoms = parse(Int32,filecontent[1])
    mol = ""
    for i in 1:natoms
        mol = mol*filecontent[2+i]
        if(i!=natoms)
            mol = mol*"\n"
        end
    end

    close(f)

    return mol

end


function molecule(mol,basis;spherical=false)

    content = split(mol,"\n")
    firstline = content[1]
    if size(split(firstline))[1] == 2
        mol = ""
	for i in 2:size(content)[1]
            mol = mol*content[i]
            if(i!=size(content)[1])
                mol = mol*"\n"
            end
        end
        charge,mul = split(firstline)
        mul = parse(Int64, mul)
        charge = parse(Int64, charge)
    else
	mul = 1
	charge = 0
    end

    if spherical
        bset = BasisSet(basis, mol)
    else
        bset = BasisSet(basis, mol, spherical=false,lib=:acsint)
    end
    p = Param(bset,mul,charge)
    p.mol = mol
    p.spherical = spherical

    return bset,p

end
