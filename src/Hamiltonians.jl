"""
    bosonic_hamiltonian(model, q_rlu; check=true, hermitian_tol=1e-10)

Assemble the dense quadratic bosonic Hamiltonian at reciprocal-lattice-unit
wave vector `q_rlu`.
"""
function bosonic_hamiltonian(
    model::SpinModel,
    q_rlu::AbstractVector{<:Real};
    check::Bool=true,
    hermitian_tol::Real=1e-10,
)
    check && validate(model)
    q = _finite_vector3(q_rlu, "q vector")
    nsites = length(model.sites)
    H = zeros(ComplexF64, 2 * nsites, 2 * nsites)
    frames = [local_frame(site.moment) for site in model.sites]
    transverse = [_transverse(frame) for frame in frames]

    for bond in model.bonds
        site_i = model.sites[bond.source]
        site_j = model.sites[bond.target]
        J = bond.scale .* model.matrices[bond.matrix].matrix
        displacement = site_j.position .+ bond.offset .- site_i.position
        phase = cis(-2pi * dot(q, displacement))
        prefactor = sqrt(site_i.spin * site_j.spin) / 2

        ui = transverse[bond.source]
        uj = transverse[bond.target]
        zi = frames[bond.source].z
        zj = frames[bond.target].z

        h11 = prefactor * _bilinear(conj.(ui), J, uj) * phase
        h12 = prefactor * _bilinear(conj.(ui), J, conj.(uj)) * phase
        h21 = prefactor * _bilinear(ui, J, uj) * phase
        h22 = prefactor * _bilinear(ui, J, conj.(uj)) * phase

        i = bond.source
        j = bond.target
        N = nsites

        H[i, j] += h11
        H[i, j + N] += h12
        H[i + N, j] += h21
        H[i + N, j + N] += h22

        H[j, i] += conj(h11)
        H[j, i + N] += conj(h21)
        H[j + N, i] += conj(h12)
        H[j + N, i + N] += conj(h22)

        longitudinal = real(_bilinear(zi, J, zj))
        H[i, i] -= site_j.spin * longitudinal
        H[j, j] -= site_i.spin * longitudinal
        H[i + N, i + N] -= site_j.spin * longitudinal
        H[j + N, j + N] -= site_i.spin * longitudinal
    end

    residual = norm(H - H')
    if residual > hermitian_tol
        throw(ArgumentError("assembled bosonic Hamiltonian is not Hermitian; residual=$residual"))
    end
    return Hermitian(H).data
end

_bilinear(a, M, b) = sum(a .* (M * b))
