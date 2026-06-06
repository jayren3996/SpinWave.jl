"""
    SpinWaveSpectrum

Raw linear spin-wave result.

Fields use stable array dimensions:

- `energies[mode, iq]`
- `correlations[a, b, mode, iq]`
"""
struct SpinWaveSpectrum
    energies::Matrix{Float64}
    correlations::Array{ComplexF64,4}
    qpoints::QPath
    diagnostics::Vector{DiagonalizationResult}
end

"""
    EnergyGrid

Broadened intensity on an energy grid with `intensity[iω, iq]`.
"""
struct EnergyGrid
    omegas::Vector{Float64}
    qpoints::QPath
    intensity::Matrix{Float64}
end

"""
    spinwave(model, path; check=true)

Compute linear spin-wave mode energies and component-resolved correlations for
all q-points in `path`.
"""
function spinwave(model::SpinModel, path::QPath; check::Bool=true)
    check && validate(model)
    nmode = length(model.sites)
    nq = length(path)
    energies = Matrix{Float64}(undef, nmode, nq)
    correlations = zeros(ComplexF64, 3, 3, nmode, nq)
    diagnostics = Vector{DiagonalizationResult}(undef, nq)
    base_corr = _site_transverse_correlations(model)

    for iq in 1:nq
        H = bosonic_hamiltonian(model, path.q_rlu[:, iq]; check=false)
        diag = diagonalize_bosonic(H)
        energies[:, iq] .= diag.energies
        diagnostics[iq] = diag
        for mode in 1:nmode
            correlations[:, :, mode, iq] .= base_corr[:, :, mode]
        end
    end

    return SpinWaveSpectrum(energies, correlations, path, diagnostics)
end

"""
    intensity(spec; components=:trace)

Return mode weights with shape `(mode, q)`.
"""
function intensity(spec::SpinWaveSpectrum; components::Symbol=:trace)
    components == :trace || throw(ArgumentError("only components=:trace is currently supported"))
    nmode, nq = size(spec.energies)
    weights = zeros(Float64, nmode, nq)
    for iq in 1:nq, mode in 1:nmode
        weights[mode, iq] = max(0.0, real(sum(spec.correlations[a, a, mode, iq] for a in 1:3)))
    end
    return weights
end

"""
    broaden(spec, omegas; eta, components=:trace)

Lorentzian-broaden a mode spectrum onto an energy grid. The returned intensity
has shape `(length(omegas), nq)`.
"""
function broaden(spec::SpinWaveSpectrum, omegas; eta::Real, components::Symbol=:trace)
    isfinite(eta) || throw(ArgumentError("eta must be finite"))
    eta > 0 || throw(DomainError(eta, "eta must be positive"))
    ωs = Float64.(collect(omegas))
    all(isfinite, ωs) || throw(ArgumentError("omegas must be finite"))
    weights = intensity(spec; components=components)
    nω = length(ωs)
    nq = size(spec.energies, 2)
    grid = zeros(Float64, nω, nq)
    for iq in 1:nq, mode in axes(spec.energies, 1)
        energy = spec.energies[mode, iq]
        weight = weights[mode, iq]
        for iω in 1:nω
            grid[iω, iq] += weight * eta / ((ωs[iω] - energy)^2 + eta^2) / pi
        end
    end
    return EnergyGrid(ωs, spec.qpoints, grid)
end

function _site_transverse_correlations(model::SpinModel)
    nmode = length(model.sites)
    corr = zeros(ComplexF64, 3, 3, nmode)
    for (i, site) in pairs(model.sites)
        frame = local_frame(site.moment)
        projector = site.spin .* (frame.x * transpose(frame.x) + frame.y * transpose(frame.y))
        corr[:, :, i] .= ComplexF64.(projector)
    end
    return corr
end
