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
    Conatins named tuples of constants. Must contain:
    "K1", "K2", "KB", and "KW".
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

    env = _resolve_environment(; Ks, temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg, m,
                               K_method, KSO4_method, BT_method, KF_method, KNH3_method,
                               Ca_method, MyAMI_mode)

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
    Conatins named tuples of constants. Must contain:
    "K1", "K2", "KB", and "KW".
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


clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

ps = (
        BT = clean(BT),
        BOH₃ = clean(BOH₃),
        BOH₄ = clean(BOH₄),
        ST   = clean(isnothing(ST) ? Constants.calc_ST(; sal).ST : ST), # else from salinity
        FT   = clean(isnothing(FT) ? Constants.calc_FT(; sal).FT : FT), # else from salinity
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT, 
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        Mg = Mg,
        Ca = Ca,
        temp_c = temp_c,
        pres_bar = pres_bar,
        sal = sal,
        pHtot = pHtot,
        pHfree = pHfree,
        pHsws = pHsws,
        pHNBS = pHNBS
    )

    # If not provided, equilibrium constants are claculated with K calculator:
    if isnothing(Ks)

        # No `kwargs...` here: the guard above has already established it is empty.
        new_Ks = K_calculator(; temp_c, sal, pres_bar, ST=nothing, FT=nothing,
        BT=nothing, K_method=K_method, KSO4_method=KSO4_method,
        BT_method=BT_method, KF_method=KF_method, KNH3_method=KNH3_method,
        Ca_method=Ca_method, MyAMI_mode=MyAMI_mode, Ca=Ca, Mg=Mg)

        ps = merge(ps, new_Ks)

    else
        Ks = Ks isa Dict ? NamedTuple(Ks) : Ks
        ps = merge(ps, (; Ks=Ks))

    end


    # Calculate pH for all scales given an input pH value & scale:
    pH_results = calc_pH_scale(
        ps.pHtot, ps.pHfree, ps.pHsws, ps.pHNBS, ps.ST, ps.FT, ps.temp_c,
        ps.sal, ps.Ks
    )
    if !isempty(pH_results)
        ps = merge(ps, pH_results)
    end

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

    pH_results_final = calc_pH_scale(
        ps.pHtot, ps.pHfree, ps.pHsws, ps.pHNBS, ps.ST, ps.FT, ps.temp_c,
        ps.sal, ps.Ks
    )
    if !isempty(pH_results_final)
        ps = merge(ps, pH_results_final)
    end

    # If any of the following parameters are known, recalculates boron isotope speciation
    has_isotope_data = !isnothing(get(ps, :ABT, nothing)) || 
                       !isnothing(get(ps, :ABOH₃, nothing)) || 
                       !isnothing(get(ps, :ABOH₄, nothing)) ||
                       !isnothing(get(ps, :δBOH₃, nothing)) || 
                       !isnothing(get(ps, :δBOH₄, nothing))
                       if has_isotope_data || !isnothing(get(ps, :δBT, nothing))
        isotope_results_final = boron_isotopes(; ps...)
        if !isempty(isotope_results_final)
            ps = merge(ps, isotope_results_final)
        end 
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
    Conatins named tuples of constants. Must contain:
    "K1", "K2", "KB", and "KW".
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

clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

ps = (
        ST   = clean(isnothing(ST) ? Constants.calc_ST(; sal).ST : ST), # else from salinity
        FT   = clean(isnothing(FT) ? Constants.calc_FT(; sal).FT : FT), # else from salinity
        BT = clean(BT),
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT, 
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        Mg = Mg,
        Ca = Ca,
        temp_c = temp_c,
        pres_bar = pres_bar,
        sal = sal,
        pHtot = pHtot,
        pHfree = pHfree,
        pHsws = pHsws,
        pHNBS = pHNBS,
        kwargs...
    )

# If not provided, equilibrium constants are claculated with K calculator:
    if isnothing(Ks)

        new_Ks = K_calculator(; temp_c, sal, pres_bar, ST=nothing, FT=nothing,
        BT=nothing, K_method=K_method, KSO4_method=KSO4_method,
        BT_method=BT_method, KF_method=KF_method, KNH3_method=KNH3_method,
        Ca_method=Ca_method, MyAMI_mode=MyAMI_mode, Ca=Ca, Mg=Mg,
        kwargs...)

        ps = merge(ps, new_Ks)

    else
        Ks = Ks isa Dict ? NamedTuple(Ks) : Ks
        ps = merge(ps, (; Ks=Ks))

    end

    # Calculate pH for all scales given an input pH value & scale:
    pH_results = calc_pH_scale(
        ps.pHtot, ps.pHfree, ps.pHsws, ps.pHNBS, ps.ST, ps.FT, ps.temp_c,
        ps.sal, ps.Ks
    )
    if !isempty(pH_results)
        ps = merge(ps, pH_results)
    end

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
    Conatins named tuples of constants. Must contain:
    "K1", "K2", "KB", and "KW".
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

    env = _resolve_environment(; Ks, temp_c, sal, pres_bar, ST, FT, BT, Ca, Mg, m,
                               K_method, KSO4_method, BT_method, KF_method, KNH3_method,
                               Ca_method, MyAMI_mode)

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

    C_count = count(!isnothing, (ps.DIC, ps.CO₂, ps.HCO₃, ps.CO₃))
    B_count = count(!isnothing, (ps.BT, ps.BOH₃, ps.BOH₄))
    iso_count = count(!isnothing, (ps.ABT,)) + count(!isnothing, (ps.ABOH₃, ps.ABOH₄))

    if !isnothing(ps.pHtot) || B_count == 2
        ps = merge(ps, B_calculator(; ps...))
        ps = merge(ps, C_calculator(; ps...))
        ps = merge(ps, calc_B_isotopes(; ps...))
    elseif iso_count == 2
        ps = merge(ps, calc_B_isotopes(; ps...))
        ps = merge(ps, B_calculator(; ps...))
        ps = merge(ps, C_calculator(; ps...))
    elseif (C_count == 2) || ((C_count == 1) && (count(!isnothing, (ps.TA, ps.BT)) == 2))
        ps = merge(ps, C_calculator(; ps...))
        ps = merge(ps, B_calculator(; ps...))
        ps = merge(ps, calc_B_isotopes(; ps...))
    else
        throw(ArgumentError("""Impossible! You haven't provided enough information.
                If you don't know pH, you must provide either:
                - Two of [DIC, CO2, HCO3, CO3] and BT
                - One of [DIC, CO2, HCO3, CO3], and TA and BT
                - Two of [BT, BO3, BO4] and one of [DIC, CO2, HCO3, CO3]
                - Two of [dBT, dBO3, dBO4] and one of [DIC, CO2, HCO3, CO3]"""))
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