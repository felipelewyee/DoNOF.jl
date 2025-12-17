const BOHR_TO_ANGSTROM = 0.5291772109
const ANGSTROM_TO_BOHR = 1.8897261246

function freq(mol::String, bset_ref, p_ref)

    # Move coordinates so center of mass is at origin
    mt, cm = center_of_mass(bset_ref)
    spin_mult, symbols, coords = string_to_xyz(mol, cm)
    mol_at_cm = xyz_data_to_string(spin_mult, symbols, coords)

    # Reference Single Point
    bset, _ = DoNOF.molecule(mol_at_cm, bset_ref.name, spherical = true)
    p = deepcopy(p_ref)
    E, C_ref, n_ref = DoNOF.energy(bset, p)

    # Numerical hessian (coords in bohr)
    coords = coords .* ANGSTROM_TO_BOHR
    hess_fin = FiniteDiff.finite_difference_jacobian(
        c -> g_linearized(c, symbols, spin_mult, C_ref, n_ref, bset, p_ref),
        coords,
    )
    #hess_fin = FiniteDiff.finite_difference_hessian(
    #    c -> e_linearized(c, symbols, spin_mult, C_ref, n_ref, bset, p_ref),
    #    coords,
    #)    

    # Mass-weighted hessian (h_ij/sqrt(m_iA m_jA))
    masses = [bset.atoms[i].mass for i = 1:bset.natoms for _ = 1:3]
    inv_sqrt_m = 1.0 ./ sqrt.(masses)
    hess_mass = (inv_sqrt_m * inv_sqrt_m') .* hess_fin

    # Explicit traslational and rotational modes
    U = zeros(3 * bset.natoms, 6)
    for iA = 1:bset.natoms
        m_root = sqrt(bset.atoms[iA].mass / mt)
        for i = 1:3
            U[3*(iA-1)+i, i] = m_root
        end
    end

    Ixyz = zeros(3, 3)
    for iA = 1:bset.natoms
        m = bset.atoms[iA].mass
        x, y, z = bset.atoms[iA].xyz
        Ixyz[1, 1] += m * (y^2 + z^2)
        Ixyz[2, 2] += m * (x^2 + z^2)
        Ixyz[3, 3] += m * (x^2 + y^2)
        Ixyz[1, 2] -= m * x * y
        Ixyz[1, 3] -= m * x * z
        Ixyz[2, 3] -= m * y * z
    end
    Ixyz = Symmetric(Ixyz)

    val, vec = eigen(Ixyz)
    for i = 1:3
        denom = sqrt(max(val[i], 1e-10))
        for iA = 1:bset.natoms
            idx = (3*(iA-1)+1):(3*iA)
            r = bset.atoms[iA].xyz
            U[idx, i+3] = sqrt(bset.atoms[iA].mass) * cross(r, vec[:, i]) / denom
        end
    end

    # Projection operator
    Id = I(3 * bset.natoms)
    P = Id - U * U'
    hess_purified = Symmetric(P * hess_mass * P)

    freq_to_molden(symbols, coords, bset, Symmetric(hess_purified), p_ref.title*".molden")

    return Symmetric(hess_purified)

end

function e_linearized(x, symbols, spin_mult, C, n, bset_ref, p_ref)

    x = x * BOHR_TO_ANGSTROM

    mol = xyz_data_to_string(spin_mult, symbols, x)

    bset, _ = DoNOF.molecule(mol, bset_ref.name, spherical = true)

    p = deepcopy(p_ref)

    E, C, n = DoNOF.energy(bset, p, C = C, n = n)

    return E
end

function g_linearized(x, symbols, spin_mult, C, n, bset_ref, p_ref)

    x = x * BOHR_TO_ANGSTROM

    mol = xyz_data_to_string(spin_mult, symbols, x)

    bset, _ = DoNOF.molecule(mol, bset_ref.name, spherical = true)

    p = deepcopy(p_ref)

    g, E, C, n = DoNOF.energy(bset, p, C = C, n = n, do_gradients = true)

    return vec(g')
end

function center_of_mass(bset)

    masses = [a.mass for a in bset.atoms]
    mt = sum(masses)
    cm = sum(bset.atoms[i].mass .* bset.atoms[i].xyz for i = 1:bset.natoms) ./ mt
    return mt, cm

end

function xyz_data_to_string(spin_mult, symbols, coords)

    natoms = length(symbols)
    io = IOBuffer()
    println(io, spin_mult)

    for iA = 1:natoms
        @printf(
            io,
            "%s % .18f % .18f % .18f",
            symbols[iA],
            coords[3*(iA-1)+1],
            coords[3*(iA-1)+2],
            coords[3*(iA-1)+3]
        )
        println(io)
    end

    return String(take!(io))
end

function string_to_xyz(mol::String, xyz)

    lines = split(strip(mol), "\n")
    lines = filter(!isempty, lines)

    spin_mult = lines[1]

    symbols = String[]
    coords = Float64[] # cm displaced to (0,0,0)

    for line in lines[2:end]
        (element, x, y, z) = split(line)
        x_val, y_val, z_val = parse.(Float64, (x, y, z))
        push!(symbols, element)
        push!(coords, x_val - xyz[1], y_val - xyz[2], z_val - xyz[3])
    end

    return spin_mult, symbols, coords
end

function freq_to_molden(symbols, coords, bset, hess, filename)

    vals, vecs = eigen(hess)

    println("freq^2")
    println(vals)
    #vals[vals .< 1e-2] .= 0
    vals[abs.(vals) .< 1e-2] .= 0

    vals = sqrt.(vals) .* 2194.746

    open(filename, "w") do io

        println(io, "[Molden Format]")
        println(io, "")
        println(io, "[ATOMS] (AU)")
        for iA = 1:bset.natoms
            @printf(
                io,
                "%s %d %d        % .10f  % .10f  % .10f",
                symbols[iA],
                iA,
                bset.atoms[iA].Z,
                coords[3*(iA-1)+1],
                coords[3*(iA-1)+2],
                coords[3*(iA-1)+3]
            )
            println(io)
        end
        println(io, "[END]")
        println(io, "[FREQ]")
        for v in vals
            #if(v > -500)
            @printf(io, "        % .10f", v)
            println(io)
            #end
        end

        println(io, "[FR-COORD]")
        for iA = 1:bset.natoms
            @printf(
                io,
                "%s          % .10f % .10f % .10f",
                symbols[iA],
                coords[3*(iA-1)+1],
                coords[3*(iA-1)+2],
                coords[3*(iA-1)+3]
            )
            println(io)
        end

        println(io, "[FR-NORM-COORD]")
        j = 0
        for (i, v) in enumerate(vals)
            #if(v > -500)
            j = j + 1
            println(io, "vibration $j")
            for iA = 1:bset.natoms
                @printf(
                    io,
                    "          % .10f % .10f % .10f",
                    vecs[3*(iA-1)+1, i],
                    vecs[3*(iA-1)+2, i],
                    vecs[3*(iA-1)+3, i]
                )
                println(io)
            end
            #end
        end

    end

end
