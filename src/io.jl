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
    Fermi.Options.set("molstring", mol)
end

function molecule_from_file(filename)
    mol = DoNOF.read_xyz_from_file(filename)

    Fermi.Options.set("molstring", mol)
end

