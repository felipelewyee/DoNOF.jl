const BOHR_TO_ANGSTROM = 0.5291772109
const ANGSTROM_TO_BOHR = 1.8897261246

function optgeo(mol::String, bset_ref, p_ref)

    # Move coordinates so center of mass is at origin

    # Reference Single Point

    # Optimize geometry (coords in bohr)
    coords = coords .* ANGSTROM_TO_BOHR
    ###### INSERTAR OPTIMIZACION ######
    # ENTRA: COORDS
    # OPTIMIZE
    # ENERGIA
        #c -> g_linearized(c, symbols, spin_mult, C_ref, n_ref, bset, p_ref)
    # GRADIENTE
        #c -> e_linearized(c, symbols, spin_mult, C_ref, n_ref, bset, p_ref)
    # SALE: RES.minimizer
    
    # coords de Bohr a Angstrom
    coords = res.minimizer * BOHR_TO_ANGSTROM
    coords[abs.(coords) .< 10^-4] .= 0

    # convertir datos a string de xyz
    # ej.
    # spin_mult: "0 1"
    # symbols: ["H", "H"]
    # coords: [x1, y1, z1, x2, y2, z2]
    mol_optimized = xyz_data_to_string(spin_mult, symbols, coords)

    return mol_optimized

end
