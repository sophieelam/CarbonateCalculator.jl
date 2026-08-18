module CarbonateCalculator
# Packages & helpers:
include("carbon.jl")
using .Carbon
include("boron.jl")
using .Boron
include("boron_isotopes.jl")
using .Isotopes
include("helpers.jl")
using .Helpers
include("constants.jl")
using .Constants
using Printf
using ForwardDiff, LinearAlgebra, Roots
include("errors.jl")
include("result.jl")
include("core.jl")
include("validation.jl")
include("display.jl")
include("conditions.jl")


export whole_system, carbon_system, boron_system, boron_isotopes, carbon_boron_calculator, carbon_calculator, propagate_errors, recalculate_at_target_conditions # export user-facing functions

# Carbon species calculations
"""
Calculates the carbon chemistry of seawater from given input parameters for 
specified output conditions.
FOR EXAMPLE: Input conditions would be the values recorded when a sample is 
measured (i.e., at sea surface level and 22C) which would be calculated for 
out put conditions representative of where the sample was taken from (i.e., 
at 500m depth and 10C).
Constants calculated using Kgen (Hain, et al., 2015)
Speciation calculations from Zeebe & Wolf-Gladrow (2001, Appendix B)

Concentration Units
-------------------
* Ca²⁺ and Mg²⁺ must be given in molar units.
* All other units must be the same and can be specified in the "unit" variable.

Parameters
----------
* pH, DIC, CO₂, HCO₃, CO₃, TA, ΩA, ΩC: array-like
    Carbon system parameters. Two must be provided to calculate the remaining.
* BT: array-like
    Total boron at the input salinity (in μ/kg). Used in total alkalinity 
    calculations. If missing, calculated from salinity: 0.0004157 * S/35.
    (Uppstrom et al. 1974)
* Ca, Mg: array-like
    The [Ca²⁺] and [Mg²⁺] of standard seawater (i.e. 35 salinity), in mol/kg. Used
    to correct the speciation constants for non-modern seawater composition via
    MyAMI, and — for Ca — to calculate saturation state.
    Both default to `nothing`, meaning modern seawater: the same composition Kgen
    assumes for the MyAMI correction, so the Ca used for saturation state and the Ca
    used to correct the constants agree. Set Ca_method to "Culkin" or "RT67" for a
    salinity-scaled concentration instead; an explicit Ca overrides Ca_method entirely.
* MyAMI_mode: str
    How the MyAMI composition correction is evaluated. Only "approximate", the
    polynomial approximation, is available; "calculate" (the full MyAMI model) is
    implemented only in the Python version of Kgen and raises an error here.
* temp_c, sal: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating speciation constants.
* pres_bar: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 

    To calculate the system at the conditions the sample was *collected* at, rather than
    those it was measured at, pass this function's result to
    `recalculate_at_target_conditions`.
* units: str
    Concentration units for carbon and boron species passed by user. All must be
    in the same unit. Can be:
    "mol", "mmol", "umol", "nmol", "pmol", or "fmol".
    Default is "umol".
* Ks: NamedTuple
    A complete seawater environment, as returned by `K_calculator`: the equilibrium
    constants in a `Ks` field, together with the totals (ST, FT, BT, Ca, Mg), fH and
    the pH-scale factors computed alongside them. Supplying it skips the constant
    calculation entirely.

        env = K_calculator(; temp_c = 25.0, sal = 35.0, pres_bar = 0.0)
        carbon_system(; TA = 2300.0, DIC = 2000.0, Ks = env)

    Not the bare constants bundle. The constants and the totals must describe the
    same water, because the pH-scale factors are built from both, so `env.Ks` alone
    is rejected rather than silently combined with salinity-derived totals. This
    entry previously said "must contain K1, K2, KB and KW", which is what led all
    four field-data scripts to pass `env.Ks` and fail on every row.
    If none, Ks are calculated with the MyAMI model. Alternative Ks for non-
    seawater conditions are available in predefined NamedTuples. See file 
    "Constants" for details.
Supplying parameters in bulk
----------------------------
Julia does this natively, so there is no `pdict` argument. Splat a NamedTuple or a Dict of
`Symbol` keys, and later entries win:

    settings = (temp_c = 2.0, sal = 34.5, K_method = "Lueker 2000")
    carbon_system(; TA = 2300.0, DIC = 2000.0, settings...)

Returns
-------
NamedTuple containing all calculated parameters
"""
function carbon_system_core(;
    pH=nothing, pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing,
    CO₃=nothing, pCO₂=nothing, fCO₂=nothing, BT=nothing, Ca=nothing,
    Mg=nothing, temp_c=25.0, sal=35.0,
    pres_bar=0.0, PT=0.0, SiT=0.0, H2ST=0.0, NH4T=0.0, ST=nothing, FT=nothing,
    pHsws=nothing, pHfree=nothing, pHNBS=nothing, unit="umol", scale="total", Ks=nothing,
    ΩC=nothing, ΩA=nothing, MyAMI_mode="approximate",
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default",
    Ca_method="default", kwargs...
)
    m = _unit_multiplier(unit)

    # BT is converted here with every other concentration, not inside the environment: it is
    # input sanitisation, and it also matches K_calculator's own convention, whose BT_method
    # fallbacks all return mol/kg.
    BT = _clean(_to_mol(BT, m))

    # K_calculator returns the complete environment - constants, totals, and the pH scale
    # factors beside them - so a supplied `Ks` is used as-is and nothing is recombined.
    env = isnothing(Ks) ? K_calculator(; temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg,
                                       K_method, KSO4_method, BT_method, KF_method,
                                       KNH3_method, Ca_method, MyAMI_mode) : Ks

    pHtot = _to_total_scale(pH, pHtot, pHsws, pHfree, pHNBS, scale, env)

    ps = (;
        kwargs...,
        DIC = _clean(_to_mol(DIC, m)),
        TA = _to_mol(TA, m),
        CO₂ = _clean(_to_mol(CO₂, m)),
        HCO₃ = _clean(_to_mol(HCO₃, m)),
        CO₃ = _clean(_to_mol(CO₃, m)),
        PT = _clean(_to_mol(PT, m)),
        SiT = _clean(_to_mol(SiT, m)),
        pCO₂ = _clean(_to_atm(pCO₂)),
        fCO₂ = _clean(_to_atm(fCO₂)),
        BT = env.BT,
        ST = env.ST,
        FT = env.FT,
        Mg = env.Mg,
        Ca = env.Ca,
        H2ST = _clean(_to_mol(H2ST, m)),
        NH4T = _clean(_to_mol(NH4T, m)),
        temp_c = temp_c,
        pres_bar = pres_bar,
        sal = sal,
        pH = nothing,     # Explicitly clear generic pH
        pHtot = pHtot,
        pHfree = nothing, # Cleared to prevent inner solvers from touching them
        pHsws = nothing,
        pHNBS = nothing,
        scale = "total",  # Force the internal scale to be total
        unit = unit,
        Ks = env.Ks
    )

    ps = merge(ps, _resolve_gases(ps))
    ps = merge(ps, _omega_to_CO₃(ps, ΩA, ΩC))

    ps = merge(ps, C_calculator(; ps...))

    ps = merge(ps, _from_total_scale(ps.pHtot, env))

    rf = calc_revelle_factor(ps.TA, ps.DIC, ps.BT, ps.PT, ps.SiT, ps.ST,
                             ps.FT, ps.H2ST, ps.NH4T, ps.Ks)
    ps = merge(ps, (revelle_factor = rf,))
    ps = merge(ps, _saturation_states(ps))

    ps = _rescale_to_unit(ps, m)
    ps = merge(ps, _rescale_gases(ps))

    return _finalise(ps)
end


# Boron species calculations
"""
Calculates the boron chemistry of seawater from given parameters. 
Constants calculated using Kgen (Hain, et al., 2015)
Speciation calculations from CBsyst (Branson, 2017)

Concentration Units
-------------------
* Ca²⁺ and Mg²⁺ must be given in molar units.
* All other units must be the same across species.

Parameters
----------
* pH, BT, BOH₃, BOH₄: array-like
    Boron system parameters. Two must be provided to calculate the remaining.
* ABT, ABOH₃, ABOH₄, δBT, δBOH₃, δBOH₄: array-like
    delta (δ) or fractional abundance (A) values for the Boron isotope system.
    One of these must be provided.
* alphaB: array-like
    The alpha value for B(OH)₃ and B(OH)₄ isotope fractionation. Default is 1.0272.
    (Kolchko, et al., 2006)
* Ca, Mg: array-like
    The [Ca²⁺] and [Mg²⁺] of standard seawater (i.e. 35 salinity), in mol/kg. Used
    to correct the speciation constants for non-modern seawater composition via
    MyAMI, and — for Ca — to calculate saturation state.
    Both default to `nothing`, meaning modern seawater: the same composition Kgen
    assumes for the MyAMI correction, so the Ca used for saturation state and the Ca
    used to correct the constants agree. Set Ca_method to "Culkin" or "RT67" for a
    salinity-scaled concentration instead; an explicit Ca overrides Ca_method entirely.
* MyAMI_mode: str
    How the MyAMI composition correction is evaluated. Only "approximate", the
    polynomial approximation, is available; "calculate" (the full MyAMI model) is
    implemented only in the Python version of Kgen and raises an error here.
* temp_c, sal: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* pres_bar: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 
* Ks: NamedTuple
    A complete seawater environment, as returned by `K_calculator`: the equilibrium
    constants in a `Ks` field, together with the totals (ST, FT, BT, Ca, Mg), fH and
    the pH-scale factors computed alongside them. Supplying it skips the constant
    calculation entirely.

        env = K_calculator(; temp_c = 25.0, sal = 35.0, pres_bar = 0.0)
        carbon_system(; TA = 2300.0, DIC = 2000.0, Ks = env)

    Not the bare constants bundle. The constants and the totals must describe the
    same water, because the pH-scale factors are built from both, so `env.Ks` alone
    is rejected rather than silently combined with salinity-derived totals. This
    entry previously said "must contain K1, K2, KB and KW", which is what led all
    four field-data scripts to pass `env.Ks` and fail on every row.
    If none, Ks are calculated with teh MyAMI model. Alternative Ks for non-
    seawater conditions are available in predefined NamedTuples. See file 
    "Constants" for details.
Supplying parameters in bulk
----------------------------
Julia does this natively, so there is no `pdict` argument. Splat a NamedTuple or a Dict of
`Symbol` keys, and later entries win:

    settings = (temp_c = 2.0, sal = 34.5, K_method = "Lueker 2000")
    carbon_system(; TA = 2300.0, DIC = 2000.0, settings...)

Returns
-------
NamedTuple containing all calculated parameters
"""
function boron_system(;
    pHtot=nothing, BT=nothing, BOH₃=nothing, BOH₄=nothing, ABT=nothing,
    ABOH₃=nothing, ABOH₄=nothing, δBT=nothing, δBOH₃=nothing, δBOH₄=nothing,
    alphaB=nothing, temp_c=25.0, sal=35.0, pres_bar =0.0, Ca=nothing, Mg=nothing, 
    ST=nothing, FT=nothing, pHsws=nothing, pHfree=nothing, 
    pHNBS=nothing, Ks=nothing, MyAMI_mode="approximate",
    K_method="default", KSO4_method="default", BT_method="default",
    KF_method="default", KNH3_method="default", Ca_method="default", kwargs...)

_reject_unknown_arguments(kwargs, boron_system)

# Check for adequate parameter input from user:
if isnothing(BT) && isnothing(BOH₃) && isnothing(BOH₄)
    throw(ArgumentError("""One of the following must be provided:
    BT, BOH₃, BOH₄"""))
end

# Check for adequate parameter input from user: 
if isnothing(δBT) && isnothing(δBOH₃) && isnothing(δBOH₄) && isnothing(ABT) && 
    isnothing(ABOH₃) && isnothing(ABOH₄)
    throw(ArgumentError("""One of the following must be provided:
    δBT, δBOH₃, δBOH₄, ABT, ABOH₃, ABOH₄"""))
end


    # Shared with the two cores, so the constants, the totals and the pH-scale conversion
    # are resolved exactly one way in this package. This function used to compute its own
    # ST/FT from salinity, call K_calculator itself, and convert pH scales through
    # `Helpers.calc_pH_scale` - a second implementation of the same conversion, and the one
    # the pHNBS sign error lived in.
    #
    # m = 1: boron_system has no `unit` argument. Speciation is linear in BT, so BOH₃ and
    # BOH₄ come back in whatever unit BT went in.
    # No `unit` argument here, so BT is used exactly as supplied and BOH₃/BOH₄ come back in
    # the same unit - valid because the speciation is linear in BT.
    BT = _clean(BT)

    # K_calculator returns the complete environment - constants, totals, and the pH scale
    # factors beside them - so a supplied `Ks` is used as-is and nothing is recombined.
    env = isnothing(Ks) ? K_calculator(; temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg,
                                       K_method, KSO4_method, BT_method, KF_method,
                                       KNH3_method, Ca_method, MyAMI_mode) : Ks

    # A supplied BT is used as given. Without one, the salinity-derived value is reported in
    # µmol/kg, matching what whole_system defaults to.
    #
    # This is the fix for a defect that made the function ignore its own primary input:
    # `ps = merge(ps, new_Ks)` used to overwrite BT with K_calculator's salinity value, so
    # BT=415.7, BT=300.0 and BT=500.0 all returned the same answer - in a function that
    # requires one of BT/BOH₃/BOH₄ to be given.
    BT_value = isnothing(BT) ? env.BT * 1e6 : env.BT

    pHtot = _to_total_scale(nothing, pHtot, pHsws, pHfree, pHNBS, "total", env)

    ps = (
        BT = BT_value,
        BOH₃ = _clean(BOH₃),
        BOH₄ = _clean(BOH₄),
        ST = env.ST,
        FT = env.FT,
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT,
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        Mg = env.Mg,
        Ca = env.Ca,
        fH = env.fH,
        temp_c = temp_c,
        pres_bar = pres_bar,
        sal = sal,
        pHtot = pHtot,
        pHfree = nothing,
        pHsws = nothing,
        pHNBS = nothing,
        Ks = env.Ks
    )

   # If pH is unknown, assign δBT value to calculate ABT, ABOH₃, and ABOH₄
    if isnothing(ps.pHtot)
        δBT = get(ps, :δBT, Isotopes.get_δBT())
        alphaB = get(ps, :alphaB, Isotopes.get_alphaB())

        # If ABT is unknown, calculate from δBT
        ABT   = !isnothing(get(ps, :ABT, nothing))   ? ps.ABT   : Isotopes.δ11_to_A11(δBT)
        
        # If δBOH₃ is known, calculate ABOH₃
        val_δBOH₃ = get(ps, :δBOH₃, nothing)
        ABOH₃ = !isnothing(val_δBOH₃) ? Isotopes.δ11_to_A11(val_δBOH₃) : nothing
        
        # If δBOH₄ is known, calculate ABOH₄
        val_δBOH₄ = get(ps, :δBOH₄, nothing)
        ABOH₄ = !isnothing(val_δBOH₄) ? Isotopes.δ11_to_A11(val_δBOH₄) : nothing

        ps = merge(ps, (; δBT, alphaB, ABT, ABOH₃, ABOH₄))

        # Calculate boron speciation and isotopes as well as pHtot
        isotope_results = calc_B_isotopes(; ps...)
        ps = merge(ps, isotope_results)
    end

    species_results = B_calculator(; ps...)
    ps = merge(ps, species_results)

    ps = merge(ps, _from_total_scale(ps.pHtot, env))

    # Solve the isotope system and report it in both notations.
    #
    # This used to re-enter the public `boron_isotopes(; ps...)`, which resolved the whole
    # environment a second time and handed `ps` - whose `Ks` is the bare constants bundle -
    # to a function whose `Ks` means the environment. Done here instead: δ in, solve, δ out.
    has_isotope_data = any(!isnothing, (ps.ABT, ps.ABOH₃, ps.ABOH₄, ps.δBOH₃, ps.δBOH₄))

    if has_isotope_data || !isnothing(ps.δBT)
        as_abundance(δ, A) = !isnothing(A) ? A :
                             (isnothing(δ) ? nothing : Isotopes.δ11_to_A11(δ))

        ps = merge(ps, (;
            alphaB = something(ps.alphaB, Isotopes.get_alphaB()),
            ABT = as_abundance(ps.δBT, ps.ABT),
            ABOH₃ = as_abundance(ps.δBOH₃, ps.ABOH₃),
            ABOH₄ = as_abundance(ps.δBOH₄, ps.ABOH₄),
        ))
        ps = merge(ps, calc_B_isotopes(; ps...))

        ps = merge(ps, (;
            δBT = !isnothing(δBT) ? δBT : A11_to_δ11(ps.ABT),
            δBOH₃ = !isnothing(δBOH₃) ? δBOH₃ : A11_to_δ11(ps.ABOH₃),
            δBOH₄ = !isnothing(δBOH₄) ? δBOH₄ : A11_to_δ11(ps.ABOH₄),
        ))
    end


    return ps

end


# Boron isotopes calculations
"""
Calculates the boron isotope chemistry of seawater from given parameters. 
Constants calculated using Kgen (Hain, et al., 2015)
Speciation calculations from CBsyst (Branson, 2017)

Concentration Units
-------------------
* Ca²⁺ and Mg²⁺ must be given in molar units.
* All other units must be the same across species.

Parameters
----------
* pH, BT, BOH₃, BOH₄: array-like
    Boron system parameters. Two must be provided to calculate the remaining.
* ABT, ABOH₃, ABOH₄, δBT, δBOH₃, δBOH₄: array-like
    delta (δ) or fractional abundance (A) values for the Boron isotope system.
    One of these must be provided.
* alphaB: array-like
    The alpha value for B(OH)₃ and B(OH)₄ isotope fractionation. Default is 1.0272.
    (Kolchko, et al., 2006)
* Ca, Mg: array-like
    The [Ca²⁺] and [Mg²⁺] of standard seawater (i.e. 35 salinity), in mol/kg. Used
    to correct the speciation constants for non-modern seawater composition via
    MyAMI, and — for Ca — to calculate saturation state.
    Both default to `nothing`, meaning modern seawater: the same composition Kgen
    assumes for the MyAMI correction, so the Ca used for saturation state and the Ca
    used to correct the constants agree. Set Ca_method to "Culkin" or "RT67" for a
    salinity-scaled concentration instead; an explicit Ca overrides Ca_method entirely.
* MyAMI_mode: str
    How the MyAMI composition correction is evaluated. Only "approximate", the
    polynomial approximation, is available; "calculate" (the full MyAMI model) is
    implemented only in the Python version of Kgen and raises an error here.
* temp_c, sal: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* pres_bar: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 
* Ks: NamedTuple
    A complete seawater environment, as returned by `K_calculator`: the equilibrium
    constants in a `Ks` field, together with the totals (ST, FT, BT, Ca, Mg), fH and
    the pH-scale factors computed alongside them. Supplying it skips the constant
    calculation entirely.

        env = K_calculator(; temp_c = 25.0, sal = 35.0, pres_bar = 0.0)
        carbon_system(; TA = 2300.0, DIC = 2000.0, Ks = env)

    Not the bare constants bundle. The constants and the totals must describe the
    same water, because the pH-scale factors are built from both, so `env.Ks` alone
    is rejected rather than silently combined with salinity-derived totals. This
    entry previously said "must contain K1, K2, KB and KW", which is what led all
    four field-data scripts to pass `env.Ks` and fail on every row.
    If none, Ks are calculated with teh MyAMI model. Alternative Ks for non-
    seawater conditions are available in predefined NamedTuples. See file 
    "Constants" for details.
Supplying parameters in bulk
----------------------------
Julia does this natively, so there is no `pdict` argument. Splat a NamedTuple or a Dict of
`Symbol` keys, and later entries win:

    settings = (temp_c = 2.0, sal = 34.5, K_method = "Lueker 2000")
    carbon_system(; TA = 2300.0, DIC = 2000.0, settings...)

Returns
-------
NamedTuple containing all calculated parameters
"""
function boron_isotopes(;
    pHtot=nothing, BT=nothing, BOH₃=nothing, BOH₄=nothing, ABT=nothing, 
    ABOH₃=nothing, ABOH₄=nothing, δBT=nothing, δBOH₃=nothing, δBOH₄=nothing, 
    alphaB=nothing, temp_c=25.0,  sal=35.0, pres_bar=0.0, Ca=nothing, Mg=nothing,
    ST=nothing, FT=nothing, pHsws=nothing, pHfree=nothing, pHNBS=nothing, 
    Ks=nothing, MyAMI_mode="approximate",
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default", Ca_method="default",
    kwargs...
)

    # Shared with the cores; see the note in boron_system. m = 1 for the same reason: no
    # `unit` argument, and the isotope ratios are dimensionless anyway.
    # No `unit` argument here, so BT is used exactly as supplied and BOH₃/BOH₄ come back in
    # the same unit - valid because the speciation is linear in BT.
    BT = _clean(BT)

    # K_calculator returns the complete environment - constants, totals, and the pH scale
    # factors beside them - so a supplied `Ks` is used as-is and nothing is recombined.
    env = isnothing(Ks) ? K_calculator(; temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg,
                                       K_method, KSO4_method, BT_method, KF_method,
                                       KNH3_method, Ca_method, MyAMI_mode) : Ks

    BT_value = isnothing(BT) ? env.BT * 1e6 : env.BT

    pHtot = _to_total_scale(nothing, pHtot, pHsws, pHfree, pHNBS, "total", env)

    ps = (
        ST = env.ST,
        FT = env.FT,
        BT = BT_value,
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT,
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        Mg = env.Mg,
        Ca = env.Ca,
        fH = env.fH,
        temp_c = temp_c,
        pres_bar = pres_bar,
        sal = sal,
        pHtot = pHtot,
        pHfree = nothing,
        pHsws = nothing,
        pHNBS = nothing,
        Ks = env.Ks,
        kwargs...
    )

    # If δBT is known, calculates ABT
    val_ABT   = !isnothing(ABT) ? ABT : (!isnothing(δBT) ? Isotopes.δ11_to_A11(δBT) : nothing)
    # If δBOH₃ is known, calculates ABOH₃
    val_ABOH₃ = !isnothing(ABOH₃) ? ABOH₃ : (!isnothing(δBOH₃) ? Isotopes.δ11_to_A11(δBOH₃) : nothing)
    # If δBOH₄ is known, calculates ABOH₄
    val_ABOH₄ = !isnothing(ABOH₄) ? ABOH₄ : (!isnothing(δBOH₄) ? Isotopes.δ11_to_A11(δBOH₄) : nothing)

    ps = merge(ps, (; 
    ABT = val_ABT, 
    ABOH₃ = val_ABOH₃, 
    ABOH₄ = val_ABOH₄
    ))

    alphaB   = !isnothing(get(ps, :alphaB, nothing))   ? ps.alphaB   : Isotopes.get_alphaB()

    b_isotopes = calc_B_isotopes(; ps...)
    ps = merge(ps, b_isotopes)

    # If unknown, calculates δBT from ABT
    final_δBT   = !isnothing(δBT) ? δBT : A11_to_δ11(ps.ABT)
    # If unknown, calculates δBOH₃ from ABOH₃
    final_δBOH₃ = !isnothing(δBOH₃) ? δBOH₃ : A11_to_δ11(ps.ABOH₃)
    # If unknown, calculates δBOH₄ from ABOH₄
    final_δBOH₄ = !isnothing(δBOH₄) ? δBOH₄ : A11_to_δ11(ps.ABOH₄)

    ps = merge(ps, (;
    δBT = final_δBT,
    δBOH₃ = final_δBOH₃,
    δBOH₄ = final_δBOH₄,
    ))

    # Report the other three pH scales. `_to_total_scale` reduced whatever the caller gave to
    # the total scale on the way in; this is the matching step on the way out.
    ps = merge(ps, _from_total_scale(ps.pHtot, env))

    return ps

end 


# Carbon & boron calculations
"""
Calculates the carbon and boron species as well as boron isotopes of seawater 
from given input parameters for specified output conditions.
FOR EXAMPLE: Input conditions would be the values recorded when a sample is 
measured (i.e., at sea surface level and 22C) which would be calculated for 
out put conditions representative of where the sample was taken from (i.e., 
at 500m depth and 10C).
Constants calculated using Kgen (Hain, et al., 2015)
Speciation calculations from Zeebe & Wolf-Gladrow (2001, Appendix B)

Note: Special Case! If pH is not known, you must provide either:
* Two of [DIC, CO₂, HCO₃, CO₃], and one of [BT, BOH₃, BOH₄]
* One of [DIC, CO₂, HCO₃, CO₃], and TA and BT
* Two of [BT, BOH₃, BOH₄] and one of [DIC, CO₂, HCO₃, CO₃]

Isotopes will only be calculated if one of [ABT, ABOH₃, ABOH₄, δBT, δBOH₃, δBOH₄]
is provided.


Concentration Units
-------------------
* Ca²⁺ and Mg²⁺ must be given in molar units.
* All other units must be the same and can be specified in the "unit" variable.
* Isotopes can be in A (11B / BT) or d (delta). Either specified, both returned.

Parameters
----------
* pH, DIC, CO₂, HCO₃, CO₃, TA, ΩA, ΩC: array-like
    Carbon system parameters. Two must be provided to calculate the remaining.
* BT: array-like
    Total boron at the input salinity (in μ/kg). Used in total alkalinity 
    calculations. If missing, calculated from salinity: 0.0004157 * S/35.
    (Uppstrom et al. 1974)
* ABT, ABOH₃, ABOH₄, δBT, δBOH₃, δBOH₄: array-like
    delta (δ) or fractional abundance (A) values for the Boron isotope system.
    One of these must be provided.
* alphaB: array-like
    The alpha value for B(OH)₃ and B(OH)₄ isotope fractionation. Default is 1.0272.
    (Kolchko, et al., 2006)
* Ca, Mg: array-like
    The [Ca²⁺] and [Mg²⁺] of standard seawater (i.e. 35 salinity), in mol/kg. Used
    to correct the speciation constants for non-modern seawater composition via
    MyAMI, and — for Ca — to calculate saturation state.
    Both default to `nothing`, meaning modern seawater: the same composition Kgen
    assumes for the MyAMI correction, so the Ca used for saturation state and the Ca
    used to correct the constants agree. Set Ca_method to "Culkin" or "RT67" for a
    salinity-scaled concentration instead; an explicit Ca overrides Ca_method entirely.
* MyAMI_mode: str
    How the MyAMI composition correction is evaluated. Only "approximate", the
    polynomial approximation, is available; "calculate" (the full MyAMI model) is
    implemented only in the Python version of Kgen and raises an error here.
* temp_c, sal: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* pres_bar: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 

    To calculate the system at the conditions the sample was *collected* at, rather than
    those it was measured at, pass this function's result to
    `recalculate_at_target_conditions`.
* units: str
    Concentration units for carbon and boron species passed by user. All must be
    in the same unit. Can be:
    "mol", "mmol", "umol", "nmol", "pmol", or "fmol".
    Default is "umol".
* Ks: NamedTuple
    A complete seawater environment, as returned by `K_calculator`: the equilibrium
    constants in a `Ks` field, together with the totals (ST, FT, BT, Ca, Mg), fH and
    the pH-scale factors computed alongside them. Supplying it skips the constant
    calculation entirely.

        env = K_calculator(; temp_c = 25.0, sal = 35.0, pres_bar = 0.0)
        carbon_system(; TA = 2300.0, DIC = 2000.0, Ks = env)

    Not the bare constants bundle. The constants and the totals must describe the
    same water, because the pH-scale factors are built from both, so `env.Ks` alone
    is rejected rather than silently combined with salinity-derived totals. This
    entry previously said "must contain K1, K2, KB and KW", which is what led all
    four field-data scripts to pass `env.Ks` and fail on every row.
    If none, Ks are calculated with teh MyAMI model. Alternative Ks for non-
    seawater conditions are available in predefined NamedTuples. See file 
    "Constants" for details.
Supplying parameters in bulk
----------------------------
Julia does this natively, so there is no `pdict` argument. Splat a NamedTuple or a Dict of
`Symbol` keys, and later entries win:

    settings = (temp_c = 2.0, sal = 34.5, K_method = "Lueker 2000")
    carbon_system(; TA = 2300.0, DIC = 2000.0, settings...)

Returns
-------
NamedTuple containing all calculated parameters

"""
function whole_system_core(;
    pH=nothing, pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing,
    CO₃=nothing, pCO₂=nothing, fCO₂=nothing, BT =nothing, BOH₃=nothing,
    BOH₄=nothing, ABT=nothing, ABOH₃=nothing, ABOH₄=nothing, δBT=nothing,
    δBOH₃=nothing, δBOH₄=nothing, alphaB=nothing, Ca=nothing,
    Mg=nothing, temp_c=25.0, sal=35.0,
    pres_bar=0.0, PT=0.0, SiT=0.0, H2ST=0.0, NH4T=0.0, ST=nothing, FT=nothing,
    pHsws=nothing, pHfree=nothing, pHNBS=nothing, unit="umol", scale="total", Ks=nothing,
    ΩC=nothing, ΩA=nothing, MyAMI_mode="approximate",
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default", Ca_method="default",
    kwargs...
)

    m = _unit_multiplier(unit)

    # BT is converted here with every other concentration, not inside the environment: it is
    # input sanitisation, and it also matches K_calculator's own convention, whose BT_method
    # fallbacks all return mol/kg.
    BT = _clean(_to_mol(BT, m))

    # K_calculator returns the complete environment - constants, totals, and the pH scale
    # factors beside them - so a supplied `Ks` is used as-is and nothing is recombined.
    env = isnothing(Ks) ? K_calculator(; temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg,
                                       K_method, KSO4_method, BT_method, KF_method,
                                       KNH3_method, Ca_method, MyAMI_mode) : Ks

    pHtot = _to_total_scale(pH, pHtot, pHsws, pHfree, pHNBS, scale, env)

    ps = (;
        kwargs...,
        DIC = _clean(_to_mol(DIC, m)),
        TA = _to_mol(TA, m),
        CO₂ = _clean(_to_mol(CO₂, m)),
        HCO₃ = _clean(_to_mol(HCO₃, m)),
        CO₃ = _clean(_to_mol(CO₃, m)),
        PT = _clean(_to_mol(PT, m)),
        SiT = _clean(_to_mol(SiT, m)),
        pCO₂ = _clean(_to_atm(pCO₂)),
        fCO₂ = _clean(_to_atm(fCO₂)),
        BT = env.BT,
        BOH₃ = _clean(_to_mol(BOH₃, m)),
        BOH₄ = _clean(_to_mol(BOH₄, m)),
        ST = env.ST,
        FT = env.FT,
        Mg = env.Mg,
        Ca = env.Ca,
        fH = env.fH,
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT,
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        H2ST = _clean(_to_mol(H2ST, m)),
        NH4T = _clean(_to_mol(NH4T, m)),
        temp_c = temp_c,
        pres_bar = pres_bar,
        sal = sal,
        pH = nothing,     # <-- Added to explicitly clear generic pH
        pHtot = pHtot,
        pHfree = nothing, # Cleared to prevent inner solvers from touching them
        pHsws = nothing,
        pHNBS = nothing,
        scale = "total", # Force the internal scale to be total
        unit = unit,
        Ks = env.Ks
    )

    alphaB   = !isnothing(get(ps, :alphaB, nothing))   ? ps.alphaB   : Isotopes.get_alphaB()
    δBT_val = (isnothing(δBT) && isnothing(ABT)) ? Isotopes.get_δBT() : δBT
    ABT_val   = !isnothing(δBT_val)   ? Isotopes.δ11_to_A11(δBT_val)   : ABT
    ABOH₃_val = !isnothing(δBOH₃) ? Isotopes.δ11_to_A11(δBOH₃) : ABOH₃
    ABOH₄_val = !isnothing(δBOH₄) ? Isotopes.δ11_to_A11(δBOH₄) : ABOH₄

    ps = merge(ps, (δBT = δBT_val, ABT = ABT_val, ABOH₃ = ABOH₃_val, ABOH₄ = ABOH₄_val))

    ps = merge(ps, _resolve_gases(ps))
    ps = merge(ps, _omega_to_CO₃(ps, ΩA, ΩC))

    # All three calculators run whatever the input; only the order changes, because pH is
    # the unknown they share. `_first_solvable` reads the same table the input validation
    # does, so the two cannot disagree about what determines the system. The branch costs
    # nothing at run time - which parameters were supplied is type information, so the
    # counts constant-fold and this resolves during specialisation.
    solve_first = _first_solvable(ps)

    if solve_first === :boron
        ps = merge(ps, B_calculator(; ps...))
        ps = merge(ps, C_calculator(; ps...))
        ps = merge(ps, calc_B_isotopes(; ps...))
    elseif solve_first === :isotopes
        ps = merge(ps, calc_B_isotopes(; ps...))
        ps = merge(ps, B_calculator(; ps...))
        ps = merge(ps, C_calculator(; ps...))
    elseif solve_first === :carbon
        ps = merge(ps, C_calculator(; ps...))
        ps = merge(ps, B_calculator(; ps...))
        ps = merge(ps, calc_B_isotopes(; ps...))
    else
        throw(ArgumentError(_underdetermined_message()))
    end

    ps = merge(ps, _from_total_scale(ps.pHtot, env))

    final_δBT   = !isnothing(δBT) ? δBT : A11_to_δ11(ps.ABT)
    final_δBOH₃ = !isnothing(δBOH₃) ? δBOH₃ : A11_to_δ11(ps.ABOH₃)
    final_δBOH₄ = !isnothing(δBOH₄) ? δBOH₄ : A11_to_δ11(ps.ABOH₄)

    ps = merge(ps, (; δBT = final_δBT, δBOH₃ = final_δBOH₃, δBOH₄ = final_δBOH₄))

    rf = calc_revelle_factor(ps.TA, ps.DIC, ps.BT, ps.PT, ps.SiT, ps.ST,
                             ps.FT, ps.H2ST, ps.NH4T, ps.Ks)
    ps = merge(ps, (revelle_factor = rf,))
    ps = merge(ps, _saturation_states(ps))

    ps = _rescale_to_unit(ps, m)
    ps = merge(ps, _rescale_gases(ps))

    return _finalise(ps)
end




"""
    carbon_system(; errors=nothing, kwargs...)

Calculates the carbonate system parameters. 
If the `errors` NamedTuple is provided, returns propagated uncertainties.
"""
function carbon_system(;
    errors=nothing, 
    pH=nothing, pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing,
    CO₃=nothing, pCO₂=nothing, fCO₂=nothing, BT=nothing, Ca=nothing,
    Mg=nothing, temp_c=25.0, sal=35.0,
    pres_bar=0.0, PT=0.0, SiT=0.0, H2ST=0.0, NH4T=0.0, ST=nothing, FT=nothing,
    pHsws=nothing, pHfree=nothing, pHNBS=nothing, unit="umol", scale="total", Ks=nothing,
    ΩC=nothing, ΩA=nothing, MyAMI_mode="approximate",
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default",
    Ca_method="default", kwargs...
    )
    # 1. Package all the inputs into a clean NamedTuple
    _reject_unknown_arguments(kwargs, carbon_system)

    inputs_nt = (
        pH=pH, pHtot=pHtot, DIC=DIC, TA=TA, CO₂=CO₂, HCO₃=HCO₃, CO₃=CO₃,
        pCO₂=pCO₂, fCO₂=fCO₂, BT=BT, Ca=Ca, Mg=Mg, temp_c=temp_c,
        sal=sal, pres_bar=pres_bar, PT=PT, SiT=SiT,
        H2ST=H2ST, NH4T=NH4T, ST=ST, FT=FT, pHsws=pHsws, pHfree=pHfree, 
        pHNBS=pHNBS, unit=unit, scale=scale, Ks=Ks, ΩC=ΩC,
        ΩA=ΩA, MyAMI_mode=MyAMI_mode, K_method=K_method, KSO4_method=KSO4_method, 
        BT_method=BT_method, KF_method=KF_method,
        KNH3_method=KNH3_method, Ca_method=Ca_method
    )

    # Merge any overflow kwargs just to be safe
    full_inputs = merge(inputs_nt, NamedTuple(kwargs))

    _check_determinacy(full_inputs, carbon_system; require_two = true)
    _check_conditions(temp_c, sal, pres_bar)
    _check_error_names(errors, full_inputs, carbon_system)

    # Identify which of the main carbonate parameters were actually provided
    # We check the inputs_nt for non-nothing values
    master_pool = [:TA, :DIC, :pH, :pHtot, :pCO₂, :fCO₂, :CO₃, :HCO₃]
    provided = [k for k in master_pool if haskey(full_inputs, k) && !isnothing(full_inputs[k])]

    if isnothing(errors)
        res = carbon_system_core(; full_inputs...)
        return CarbonateResult(res, nothing, provided,
                               _settings(full_inputs), full_inputs, :carbon, nothing)
    else
        res_data = propagate_errors(carbon_system_core; inputs=full_inputs, errors=errors)
        return CarbonateResult(res_data.val, res_data.err, provided,
                               _settings(full_inputs), full_inputs, :carbon, errors)
    end
end


function whole_system(;
    errors=nothing, 
    pH=nothing, pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing,
    CO₃=nothing, pCO₂=nothing, fCO₂=nothing, BT=nothing, BOH₃=nothing,
    BOH₄=nothing, ABT=nothing, ABOH₃=nothing, ABOH₄=nothing, δBT=nothing,
    δBOH₃=nothing, δBOH₄=nothing, alphaB=nothing, Ca=nothing,
    Mg=nothing, temp_c=25.0, sal=35.0,
    pres_bar=0.0, PT=0.0, SiT=0.0, H2ST=0.0, NH4T=0.0, ST=nothing, FT=nothing,
    pHsws=nothing, pHfree=nothing, pHNBS=nothing, unit="umol", scale="total", Ks=nothing,
    ΩC=nothing, ΩA=nothing, MyAMI_mode="approximate",
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default", Ca_method="default",
    kwargs...
)
    
    # 1. Package all the inputs into a clean NamedTuple
    _reject_unknown_arguments(kwargs, whole_system)

    inputs_nt = (
        pH=pH, pHtot=pHtot, DIC=DIC, TA=TA, CO₂=CO₂, HCO₃=HCO₃, CO₃=CO₃,
        pCO₂=pCO₂, fCO₂=fCO₂, BT=BT, BOH₃=BOH₃, BOH₄=BOH₄, ABT=ABT,
        ABOH₃=ABOH₃, ABOH₄=ABOH₄, δBT=δBT, δBOH₃=δBOH₃, δBOH₄=δBOH₄, 
        alphaB=alphaB, Ca=Ca, Mg=Mg, temp_c=temp_c, sal=sal,
        pres_bar=pres_bar, PT=PT, SiT=SiT, H2ST=H2ST,
        NH4T=NH4T, ST=ST, FT=FT, pHsws=pHsws, pHfree=pHfree, pHNBS=pHNBS, 
        unit=unit, scale=scale, Ks=Ks, ΩC=ΩC, ΩA=ΩA,
        MyAMI_mode=MyAMI_mode, K_method=K_method, KSO4_method=KSO4_method, 
        BT_method=BT_method, KF_method=KF_method,
        KNH3_method=KNH3_method, Ca_method=Ca_method
    )

    # Merge any overflow kwargs just to be safe
    full_inputs = merge(inputs_nt, NamedTuple(kwargs))

    # require_two = false: the boron and isotope routes to a solution live in the core's own
    # check, which reports them properly.
    _check_determinacy(full_inputs, whole_system; require_two = false)
    _check_conditions(temp_c, sal, pres_bar)
    _check_error_names(errors, full_inputs, whole_system)

    # Identify which of the main carbonate parameters were actually provided
    # We check the inputs_nt for non-nothing values
    master_pool = [:TA, :DIC, :pH, :pHtot, :pCO₂, :fCO₂, :CO₃, :HCO₃]
    provided = [k for k in master_pool if haskey(full_inputs, k) && !isnothing(full_inputs[k])]

    if isnothing(errors)
        # Fast path
        res = whole_system_core(; full_inputs...)
        return CarbonateResult(res, nothing, provided,
                               _settings(full_inputs), full_inputs, :whole, nothing)
    else
        res_data = propagate_errors(whole_system_core; inputs=full_inputs, errors=errors)
        return CarbonateResult(res_data.val, res_data.err, provided,
                               _settings(full_inputs), full_inputs, :whole, errors)
    end
end


"""
    carbon_calculator(; kwargs...)

Alias for [`carbon_system`](@ref).

The two were separate implementations only because `carbon_system` also handled output
conditions; that job now belongs to [`recalculate_at_target_conditions`](@ref), so the
duplicate has nothing left to do.
"""
carbon_calculator(; kwargs...) = carbon_system(; kwargs...)

"""
    carbon_boron_calculator(; kwargs...)

Alias for [`whole_system`](@ref). See [`carbon_calculator`](@ref) for why.
"""
carbon_boron_calculator(; kwargs...) = whole_system(; kwargs...)

end # module