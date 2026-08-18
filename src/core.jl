# The pipeline shared by `carbon_system_core` and `whole_system_core`.
#
# The two cores ran the same stages in the same order and differed only in the boron and
# isotope block: 200 and 241 lines with 29 unique to the carbon one, and `carbon_system`'s
# 40 output keys are a strict subset of `whole_system`'s 50. That duplication had already
# produced two live defects by the time it was split out — `whole_system` could not take
# ΩA/ΩC as inputs, and `carbon_system` did not unit-rescale `Alk_H2S`/`Alk_NH3` — each
# present in one copy and not the other. Both are fixed here by there being one copy.
#
# The stages are named functions rather than one flag-driven core because the cores marked
# them with comments (`# --- THE WAY IN ---`, `# Converting values back into their original
# units`), which is the signal to split.

"Concentration unit names to the factor that converts them to mol/kg."
const UNIT_SCALES = Dict(
    "mol" => 1.0,
    "mmol" => 1.0e3,
    "umol" => 1.0e6,
    "nmol" => 1.0e9,
    "pmol" => 1.0e12,
    "fmol" => 1.0e15,
)

_unit_multiplier(unit) = get(UNIT_SCALES, unit, 1.0)

"A negative concentration is not a small one; it is a broken input, and NaN says so."
_clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

"Convert a concentration from the caller's unit to mol/kg."
_to_mol(v, m) = isnothing(v) ? nothing : v / m

"Convert a gas from µatm to atm."
_to_atm(v) = isnothing(v) ? nothing : v / 1e6

"""
Every concentration that is held internally in mol/kg and must be converted back to the
caller's unit on the way out.

One list for both systems. It used to be written out separately in each core, and the
carbon copy was missing `Alk_H2S` and `Alk_NH3` — so at `unit="mmol"` those two came back
in mol/kg while every other concentration in the same result was in mmol/kg, a factor of
1000 apart. Keys absent from a given system are skipped.
"""
const CONCENTRATION_KEYS = (:DIC, :TA, :BT, :CO₂, :HCO₃, :CO₃, :PT, :SiT, :BOH₃, :BOH₄,
                            :CAlk, :BAlk, :PAlk, :OH, :SiAlk, :HSO₄, :Hfree, :HF,
                            :Alk_H2S, :Alk_NH3)

"""
    _resolve_environment(; Ks, temp_c, sal, pres_bar, ...)

Resolve the equilibrium constants and the seawater composition they belong to.

Returns the constants alongside the totals and the two pH-scale conversion factors, so that
everything downstream reads one consistent description of the water rather than re-deriving
pieces of it.
"""
function _resolve_environment(; Ks, temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg, m,
                              K_method, KSO4_method, BT_method, KF_method, KNH3_method,
                              Ca_method, MyAMI_mode)
    Ks_env = if isnothing(Ks)
        K_calculator(;
            temp_c=temp_c, sal=sal, pres_bar=pres_bar, ST=ST, FT=FT, BT=BT,
            K_method=K_method, KSO4_method=KSO4_method, BT_method=BT_method,
            KF_method=KF_method, KNH3_method=KNH3_method, Ca_method=Ca_method,
            MyAMI_mode=MyAMI_mode, Ca=Ca, Mg=Mg
        )
    else
        Ks
    end

    # Ca drives both the MyAMI correction and the saturation state. K_calculator has
    # already resolved it - the explicit value if one was given, Ca_method otherwise - so
    # take it from there rather than re-deriving it and risking the two disagreeing.
    Ca_val = hasproperty(Ks_env, :Ca) ? Ks_env.Ca : something(Ca, Constants.MODERN_CALCIUM)
    Mg_val = hasproperty(Ks_env, :Mg) ? Ks_env.Mg : something(Mg, Constants.MODERN_MAGNESIUM)

    ST_val = Ks_env.ST
    FT_val = Ks_env.FT
    KS_val = Ks_env.Ks.KS
    KF_val = Ks_env.Ks.KF

    BT_from_env = hasproperty(Ks_env, :BT) ? Ks_env.BT :
                  (hasproperty(Ks_env.Ks, :BT) ? Ks_env.Ks.BT : 0.0)

    internal_BT = !isnothing(BT) ? _clean(_to_mol(BT, m)) :
                  (BT_from_env > 1.0 ? BT_from_env / 1e6 : BT_from_env)

    fH_val = Ks_env.fH isa AbstractArray ? Ks_env.fH[1] : Ks_env.fH

    tf_fac = (1.0 + ST_val / KS_val)
    ts_fac = tf_fac / (1.0 + ST_val / KS_val + FT_val / KF_val)

    return (Ks = Ks_env.Ks, ST = ST_val, FT = FT_val, BT = internal_BT,
            Ca = Ca_val, Mg = Mg_val, fH = fH_val, tf_fac = tf_fac, ts_fac = ts_fac)
end

"""
    _to_total_scale(pH, pHtot, pHsws, pHfree, pHNBS, scale, env)

Reduce whatever pH the caller gave, on whatever scale, to a single total-scale value.

Everything inside the solvers works on the total scale, so this is the one place a scale
conversion happens on the way in — and `_from_total_scale` the one place on the way out.
"""
function _to_total_scale(pH, pHtot, pHsws, pHfree, pHNBS, scale, env)
    # A generic `pH` means "on the scale named by `scale`".
    if !isnothing(pH) && isnothing(pHtot) && isnothing(pHsws) && isnothing(pHfree) &&
       isnothing(pHNBS)
        named = lowercase(scale)
        named == "total" && (pHtot = pH)
        named == "sws" && (pHsws = pH)
        named == "free" && (pHfree = pH)
        named == "nbs" && (pHNBS = pH)
    end

    if isnothing(pHtot)
        if !isnothing(pHsws)
            pHtot = pHsws - log10(env.ts_fac)
        elseif !isnothing(pHfree)
            pHtot = pHfree - log10(env.tf_fac)
        elseif !isnothing(pHNBS)
            pHtot = pHNBS + log10(env.fH) - log10(env.ts_fac)
        end
    end

    return pHtot
end

"The other three pH scales, from the solved total-scale value."
function _from_total_scale(pHtot, env)
    pHsws = pHtot + log10(env.ts_fac)
    return (pHfree = pHtot + log10(env.tf_fac),
            pHsws = pHsws,
            pHNBS = pHsws - log10(env.fH))
end

"""
Resolve CO₂, fCO₂ and pCO₂ from whichever of the three was supplied.

`fCO₂` is the one the equilibrium constant `K0` relates to CO₂; `pCO₂` reaches it through
the virial correction in `pCO₂_to_fCO₂`.
"""
function _resolve_gases(ps)
    CO₂, fCO₂, pCO₂ = ps.CO₂, ps.fCO₂, ps.pCO₂

    if isnothing(CO₂)
        if !isnothing(fCO₂)
            CO₂ = fCO₂ * ps.Ks.K0
            pCO₂ = isnothing(pCO₂) ? fCO₂_to_pCO₂(fCO₂, ps.temp_c) : pCO₂
        elseif !isnothing(pCO₂)
            fCO₂ = pCO₂_to_fCO₂(pCO₂, ps.temp_c)
            CO₂ = fCO₂ * ps.Ks.K0
        end
    end

    return (CO₂ = CO₂, fCO₂ = fCO₂, pCO₂ = pCO₂)
end

"""
Convert a supplied saturation state into the [CO₃²⁻] that produces it.

Ω is a statement about carbonate ion concentration, so this lets ΩA or ΩC stand in as one
of the two carbonate parameters. It lived only in the carbon core, which is why
`whole_system(TA=2300.0, ΩC=5.0)` used to fail with "Impossible! You haven't provided
enough information" while `carbon_system` solved the same input.
"""
function _omega_to_CO₃(ps, ΩA, ΩC)
    calcium = ps.Ca * ps.sal / 35.0
    if !isnothing(ΩA)
        return (CO₃ = ΩA * ps.Ks.KspA / calcium,)
    elseif !isnothing(ΩC)
        return (CO₃ = ΩC * ps.Ks.KspC / calcium,)
    end
    return NamedTuple()
end

"Saturation states of aragonite and calcite for the solved [CO₃²⁻]."
function _saturation_states(ps)
    calcium = ps.Ca * ps.sal / 35.0
    return (ΩA = ps.CO₃ * calcium / ps.Ks.KspA,
            ΩC = ps.CO₃ * calcium / ps.Ks.KspC)
end

"""
Convert concentrations from the internal mol/kg back to the caller's unit.

`DIC`, `TA` and `BT` are left alone when already greater than 1.0: they can be supplied
either in the caller's unit or as mol/kg, and a value above 1 mol/kg is not physical, so it
must already have been converted.
"""
function _rescale_to_unit(ps, m)
    m == 1 && return ps

    rescaled = (; [
        begin
            value = getfield(ps, k)
            if isnothing(value)
                k => nothing
            elseif (k === :BT || k === :DIC || k === :TA) && value > 1.0
                k => value
            else
                k => value * m
            end
        end for k in CONCENTRATION_KEYS if hasproperty(ps, k)
    ]...)

    return merge(ps, rescaled)
end

"Gases are held internally in atm and reported in µatm."
_rescale_gases(ps) = (pCO₂ = ps.pCO₂ * 1e6, fCO₂ = ps.fCO₂ * 1e6)

"Drop the internal bookkeeping that is not part of the result."
function _finalise(ps)
    keys_to_keep = Tuple(k for k in keys(ps) if k != :scale)
    return NamedTuple{keys_to_keep}(ps)
end
