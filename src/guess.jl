function read_C(; title = "donof")

    C = nothing
    try
        C = load(title * ".jld2")["C"]
    catch
        try
            C = load(title * ".jld")["C"]
        catch
            println(title * "_C.jld2 not found")
        end
    end
    return C

end

function read_n(; title = "donof")

    n = nothing
    try
        n = load(title * ".jld2")["n"]
    catch
        try
            n = load(title * ".jld")["n"]
        catch
            println(title * "_n.jld2 not found")
        end
    end
    return n

end

function read_fmiug0(; title = "donof")

    fmiug0 = nothing
    try
        fmiug0 = load(title * ".jld2")["fmiug0"]
    catch
        try
            fmiug0 = load(title * ".jld")["fmiug0"]
        catch
            println(title * "_fmiug0.jld2 not found")
        end
    end
    return fmiug0

end

function read_all(; title = "donof")
    C = read_C(title = title)
    n = read_n(title = title)
    fmiug0 = read_fmiug0(title = title)

    return C, n, fmiug0

end

function order_subspaces(old_C, old_n, elag, H, I, p)

    C = zeros(p.nbf, p.nbf)
    n = zeros(p.nbf5)

    #Sort no1 orbitals
    elag_diag = diag(elag)[1:p.no1]
    sort_idx = sortperm(elag_diag)
    for i = 1:p.no1
        i_old = sort_idx[i]
        C[1:p.nbf, i] = old_C[1:p.nbf, i_old]
    end

    #Sort ndoc subspaces
    elag_diag = diag(elag)[p.no1+1:p.ndoc]
    sort_idx = sortperm(elag_diag)
    C[1:p.nbf, 1:p.no1] = old_C[1:p.nbf, 1:p.no1]
    n[1:p.no1] = old_n[1:p.no1]
    for i = 1:p.ndoc
        i_old = sort_idx[i]
        ll = p.no1 + p.ndns + (p.ndoc - i) + 1
        ul = ll + p.ndoc * (p.ncwo - 1)

        ll_old = p.no1 + p.ndns + (p.ndoc - i_old) + 1
        ul_old = ll_old + p.ndoc * (p.ncwo - 1)

        C[1:p.nbf, p.no1+i] = old_C[1:p.nbf, p.no1+i_old]
        C[1:p.nbf, ll:p.ndoc:ul] = old_C[1:p.nbf, ll_old:p.ndoc:ul_old]
        n[p.no1+i] = old_n[p.no1+i_old]
        n[ll:p.ndoc:ul] = old_n[ll_old:p.ndoc:ul_old]
    end

    #Sort nsoc orbitals
    elag_diag = diag(elag)[p.no1+p.ndoc+1:p.no1+p.ndns]
    sort_idx = sortperm(elag_diag)
    for i = 1:p.nsoc
        i_old = sort_idx[i]
        C[1:p.nbf, p.no1+p.ndoc+i] = old_C[1:p.nbf, p.no1+p.ndoc+i_old]
        n[p.no1+p.ndoc+i] = old_n[p.no1+p.ndoc+i_old]
    end

    C[1:p.nbf, p.nbf5+1:p.nbf] = old_C[1:p.nbf, p.nbf5+1:p.nbf]

    cj12, ck12 = PNOFi_selector(n, p)
    Etmp, elag, sumdiff, maxdiff = ENERGY1r(C, n, H, I, cj12, ck12, p)

    return C, n, elag

end

function write_to_DoNOFsw(p, bset, n, C, elag, fmiug0, it, E)

    f = open(p.title * ".gcf", "w")
    for i = 1:p.nbf5
        @printf(f, "%6d %30.16f\n", i, n[i])
    end
    for i = p.nbf5+1:p.nbf
        @printf(f, "%6d %30.16f\n", i, 0.0)
    end
    sumsl = p.nbeta
    for i = 1:p.nbeta
        sumsl -= n[i]
    end
    @printf(f, "%30.16f\n", sumsl)
    for i = 1:p.nbf
        for j = 1:p.nbf
            @printf(f, "%6d %30.16f\n", j, C[j, i])
        end
    end
    for i = 1:p.nbf
        @printf(f, "%6d %30.16f\n", i, elag[i])
    end
    for i = 1:p.nbf
        @printf(f, "%6d %30.16f\n", i, fmiug0[i])
    end
    @printf(f, "%6d %30.16f\n", 0, E)
    @printf(f, "%6d %6d %6d %6d %6d %6d\n", p.no1, p.ndoc, p.nsoc, p.ncwo, p.nac, p.no0)
    for i = 1:bset.natoms
        @printf(
            f,
            "%6d%30.16f%30.16f%30.16f\n",
            bset.atoms[i].Z,
            bset.atoms[i].xyz[1],
            bset.atoms[i].xyz[2],
            bset.atoms[i].xyz[3]
        )
    end

    close(f)

end

function read_from_DoNOFsw(p, bset)

    f = open(p.title * ".gcf", "r")
    n = zeros(p.nbf5)
    for i = 1:p.nbf5
        line = readline(f)
        j, val = split(line)
        n[i] = parse(Float64, val)
    end
    for i = p.nbf5+1:p.nbf
        line = readline(f)
    end
    sumsl = readline(f)
    C = zeros(p.nbf, p.nbf)
    for i = 1:p.nbf
        for j = 1:p.nbf
            line = readline(f)
            k, val = split(line)
            C[j,i] = parse(Float64, val)
        end
    end
    for i = 1:p.nbf
        line = readline(f)
    end
    for i = 1:p.nbf
        line = readline(f)
    end

    close(f)

    return C,n

end

function guess_gamma_trigonometric(ndoc, ncwo)

    nv = (ncwo) * ndoc
    γ = zeros(nv)
    for i = 1:ndoc
        γ[i] = acos(sqrt(2.0 * 0.999 - 1.0))
        ll = ndoc + (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 2)
        γ_coupled = @view γ[ll:ndoc:ul]
        for j = 1:ncwo-1
            γ_coupled[j] = asin(sqrt(1.0 / (ncwo - j + 1)))
        end
    end
    return γ
end

function guess_gamma_softmax(ndoc, ncwo)

    nv = ncwo * ndoc
    γ = zeros(nv)
    for i = 1:ndoc
        ll = (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 1)
        γ_coupled = @view γ[ll:ndoc:ul]
        for j = 1:ncwo
            γ_coupled[j] = log(max(0.001 / ncwo, 1e-16)) - log(0.999)
        end
    end
    return γ
end

function guess_gamma_ebi(ndoc, nbf)
    """Compute a guess for gammas in the softmax parameterization
    of the occupation numbers"""

    #TODO: Add equations and look to reduce a variable

    nv = nbf
    γ = zeros(nv)
    for i = 1:ndoc
        γ[i] = erfinv(2 * 0.999 - 1)
    end
    val = ndoc * (1.0 - 0.999)
    for i = ndoc+1:nv
        γ[i] = erfinv(2 * val / nv - 1)
    end
    return γ
end
