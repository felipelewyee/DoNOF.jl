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

function molecule(mol)

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
        Fermi.Options.set("multiplicity", parse(Int64, mul))
        Fermi.Options.set("charge", parse(Int64, charge))
    end

    Fermi.Options.set("molstring", mol)
end

function molecule_from_file(filename)
    mol = DoNOF.read_xyz_from_file(filename)

    Fermi.Options.set("molstring", mol)
end

