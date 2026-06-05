module Calculator
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
using PythonCall
const kgen = pyimport("kgen")
using Printf
const np = pyimport("numpy")


export whole_system, carbon_system, boron_system, boron_isotopes # export user-facing functions



# Carbon species calculations
"""
Calculates the carbon chemistry of seawater from given parameters. 
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
    to calculate MyAMI constants.
* T_in, S_in: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* P_in: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 
* T_out, S_out: array-like
    Temperature in Celcius and salinity in PSU of desired output conditions. Used
    in calculating constants.
* P_out: array-like
    Pressure in bar of desired output conditions. Used in pressure-correcting
    constants.
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
* pdict: dict
    Optional: can be used to provide some or all paramters as a dictionary/
    NamedTuples with the same key names. Any paramters in pdict will 
    overwrite manually specified paramters.

Returns
-------
NamedTuple containing all calculated parameters
"""
function carbon_system(;
    pH=nothing, pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing,
    CO₃=nothing, pCO₂=nothing, fCO₂=nothing, BT=nothing, Ca=0.0102821,
    Mg=0.0528171, T_in=25.0, T_out=nothing, S_in=35.0, S_out=nothing,
    P_in=0.0, P_out=nothing, PT=0.0, SiT=0.0, H2ST=0.0, NH4T=0.0, ST=nothing, FT=nothing,
    pHsws=nothing, pHfree=nothing, pHNBS=nothing, unit="umol", scale="total", Ks=nothing,
    pdict=nothing, ΩC=nothing, ΩA=nothing, MyAMI_mode="calculate", 
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", K_mode="static", KNH3_method="default",
    Ca_method="default", kwargs...
)
    # Assigning scaling factors to units:
    udict = Dict(
        "mol" => 1.0,
        "mmol" => 1.0e3,
        "umol" => 1.0e6,
        "nmol" => 1.0e9,
        "pmol" => 1.0e12,
        "fmol" => 1.0e15
    )

    m = get(udict, unit, 1.0)
    scale_unit(v) = isnothing(v) ? nothing : v / m
    scale_gas(v) = isnothing(v) ? nothing : v / 1e6 
    clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

    # 1. Pre-calculate Ks_env using direct function arguments
    if isnothing(Ks)
        Ks_env = K_calculator(;
            T_in=T_in, S_in=S_in, P_in=P_in, ST=ST, FT=FT, BT=BT,
            K_method=K_method, KSO4_method=KSO4_method, BT_method=BT_method,
            KF_method=KF_method, KNH3_method=KNH3_method, Ca_method=Ca_method, K_mode=K_mode
        )
    else
        Ks_env = Ks
    end

    ST_val = Ks_env.ST
    FT_val = Ks_env.FT
    KS_val = Ks_env.Ks.KS
    KF_val = Ks_env.Ks.KF
    
    BT_from_env = hasproperty(Ks_env, :BT) ? Ks_env.BT : 
                  (hasproperty(Ks_env.Ks, :BT) ? Ks_env.Ks.BT : 0.0)

    internal_BT = !isnothing(BT) ? clean(scale_unit(BT)) : 
                  (BT_from_env > 1.0 ? BT_from_env / 1e6 : BT_from_env)

    fH_val = Ks_env.fH isa AbstractArray ? Ks_env.fH[1] : Ks_env.fH
    
    tf_fac = (1.0 + ST_val / KS_val) 
    ts_fac = tf_fac / (1.0 + ST_val / KS_val + FT_val / KF_val)
    
    # --- THE WAY IN: Map generic pH, then strictly convert to pHtot ---
    if !isnothing(pH) && isnothing(pHtot) && isnothing(pHsws) && isnothing(pHfree) && isnothing(pHNBS)
        if lowercase(scale) == "total"
            pHtot = pH
        elseif lowercase(scale) == "sws"
            pHsws = pH
        elseif lowercase(scale) == "free"
            pHfree = pH
        elseif lowercase(scale) == "nbs"
            pHNBS = pH
        end
    end

    if isnothing(pHtot)
        if !isnothing(pHsws)
            pHtot = pHsws - log10(ts_fac)
        elseif !isnothing(pHfree)
            pHtot = pHfree - log10(tf_fac)
        elseif !isnothing(pHNBS)
            pHtot = pHNBS + log10(fH_val) - log10(ts_fac)
        end
    end

    ps = (;
        kwargs...,
        DIC = clean(scale_unit(DIC)),
        TA = scale_unit(TA),
        CO₂ = clean(scale_unit(CO₂)),
        HCO₃ = clean(scale_unit(HCO₃)),
        CO₃ = clean(scale_unit(CO₃)),
        PT = clean(scale_unit(PT)),
        SiT = clean(scale_unit(SiT)),
        pCO₂ = clean(scale_gas(pCO₂)),
        fCO₂ = clean(scale_gas(fCO₂)),
        BT = internal_BT,
        ST = Ks_env.ST,
        FT = Ks_env.FT,
        Mg = Mg,
        Ca = Ca,
        H2ST = clean(scale_unit(H2ST)),
        NH4T = clean(scale_unit(NH4T)),
        T_in = T_in,
        P_in = P_in,
        S_in = S_in,
        T_out = T_out,
        S_out = S_out,
        P_out = P_out,
        pH = nothing,     # Explicitly clear generic pH
        pHtot = pHtot,
        pHfree = nothing, # Cleared to prevent inner solvers from touching them
        pHsws = nothing,
        pHNBS = nothing,
        scale = "total",  # Force the internal scale to be total
        unit = unit,
        Ks = Ks_env.Ks
    )

    # --- Gas Conversion Logic ---
    local_fCO2 = ps.fCO₂
    local_pCO2 = ps.pCO₂
    local_CO2  = ps.CO₂

    if isnothing(local_CO2)
        if !isnothing(local_fCO2)
            local_CO2  = local_fCO2 * ps.Ks.K0
            local_pCO2 = isnothing(local_pCO2) ? fCO₂_to_pCO₂(local_fCO2, ps.T_in) : local_pCO2
        elseif !isnothing(local_pCO2)
            local_fCO2 = pCO₂_to_fCO₂(local_pCO2, ps.T_in)
            local_CO2  = local_fCO2 * ps.Ks.K0
            local_pCO2 = local_pCO2
        end
    end

    ps = merge(ps, (CO₂ = local_CO2, fCO₂ = local_fCO2, pCO₂ = local_pCO2))

    # If ΩA or ΩC are provided, use them to calculate CO₃:
    if !isnothing(ΩA)
        new_CO₃ = ΩA * ps.Ks.KspA / (ps.Ca * ps.S_in / 35.0)
        ps = merge(ps, (CO₃ = new_CO₃,)) 
    elseif !isnothing(ΩC)
        new_CO₃ = ΩC * ps.Ks.KspC / (ps.Ca * ps.S_in / 35.0)
        ps = merge(ps, (CO₃ = new_CO₃,))
    end 

    # Calculate all of the carbon chemistry for input conditions: 
    C_calculations = C_calculator(; ps...)
    ps = merge(ps, C_calculations)

    # --- THE WAY OUT: Calculate other pH scales from the finished Total pH ---
    solved_pHtot = ps.pHtot
    final_pHsws  = solved_pHtot + log10(ts_fac)
    final_pHfree = solved_pHtot + log10(tf_fac)
    final_pHNBS  = final_pHsws - log10(fH_val) 

    ps = merge(ps, (;
        pHfree = final_pHfree,
        pHsws  = final_pHsws,
        pHNBS  = final_pHNBS
    ))

    # Calculate the Revelle Factor for input conditions:
    rf = calc_revelle_factor(ps.TA, ps.DIC, ps.BT, ps.PT, ps.SiT, ps.ST, 
                             ps.FT, ps.H2ST, ps.NH4T, ps.Ks)
    ps = merge(ps, (revelle_factor = rf,))

    # Re-calculate Ω values:
    oCa = ps.Ca * ps.S_in / 35.0
    ΩA = ps.CO₃ * oCa / ps.Ks.KspA
    ΩC = ps.CO₃ * oCa / ps.Ks.KspC
    ps = merge(ps, (ΩA = ΩA, ΩC = ΩC))

    # Converting values back into their original units:
    if m ≠ 1
        target_keys = (:DIC, :TA, :BT, :CO₂, :HCO₃, :CO₃, :PT, :SiT, :CAlk, :BAlk, :PAlk, :OH, :SiAlk, :HSO₄, :Hfree, :HF)
        
        rescaled_concs = (; [
            begin
                val = getfield(ps, k)
                if isnothing(val)
                    k => nothing
                elseif (k === :BT || k === :DIC || k === :TA) && val > 1.0
                    k => val
                else
                    k => val * m
                end
            end for k in target_keys if hasproperty(ps, k)
        ]...)
        
        ps = merge(ps, rescaled_concs)
    end

    rescaled_gases = (
        pCO₂ = ps.pCO₂ * 1e6,
        fCO₂ = ps.fCO₂ * 1e6
    )
    ps = merge(ps, rescaled_gases)

    # Assigning output conditions:
    if (!isnothing(ps.T_out)) || (!isnothing(ps.S_out)) || (!isnothing(ps.P_out))
        out_params = (
            T_out = isnothing(ps.T_out) ? ps.T_in : ps.T_out,
            S_out = isnothing(ps.S_out) ? ps.S_in : ps.S_out,
            P_out = isnothing(ps.P_out) ? ps.P_in : ps.P_out
        )
        ps = merge(ps, out_params)

        local_BT = ps.S_in ≠ ps.S_out ? ps.BT .* ps.S_out ./ ps.S_in : ps.BT
        local_ST = ps.S_in ≠ ps.S_out ? ps.ST .* ps.S_out ./ ps.S_in : ps.ST
        local_FT = ps.S_in ≠ ps.S_out ? ps.FT .* ps.S_out ./ ps.S_in : ps.FT

        # Re-calculating for different output conditions
        out_cond = carbon_system(;
            TA = ps.TA,
            DIC = ps.DIC,
            T_in = ps.T_out,
            S_in = ps.S_out,
            P_in = ps.P_out,
            T_out = nothing, 
            S_out = nothing,
            P_out = nothing,
            unit = ps.unit, #"mol", # Molar units to match internal system depth
            scale = ps.scale,
            Ca = ps.Ca,
            Mg = ps.Mg,
            BT = local_BT,
            FT = local_FT,
            ST = local_ST,
            PT = ps.PT,
            SiT = ps.SiT,
            MyAMI_mode = MyAMI_mode, 
            K_method="default",
            KSO4_method="default",
            BT_method="default"
        )

        outputs = [
            "BAlk", "BT", "CAlk", "CO₂", "CO₃", "DIC", "H", "HCO₃", "HF", "HSO₄", 
            "Hfree", "Ks", "OH", "PAlk", "SiAlk", "TA", "FT","PT", "ST", "SiT", 
            "fCO₂", "pCO₂", "pHfree", "pHsws", "pHtot", "pHNBS", "ΩA", "ΩC", 
            "revelle_factor"
        ]
        
        in_vals = (; (Symbol(k * "_in") => getfield(ps, Symbol(k)) 
                      for k in outputs if hasproperty(ps, Symbol(k)))...)

        out_vals = (; (Symbol(k) => getfield(out_cond, Symbol(k)) 
                       for k in outputs if hasproperty(out_cond, Symbol(k)))...)
        ps = merge(ps, in_vals, out_vals)
    end 

    # Clean up tuple before returning
    keys_to_keep = Tuple(k for k in keys(ps) if k ∉ (:pdict, :unit, :scale, :unit))
    ps = NamedTuple{keys_to_keep}(ps)

    return ps
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
    to calculate MyAMI constants.
* T_in, S_in: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* P_in: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 
* Ks: NamedTuple
    Conatins named tuples of constants. Must contain:
    "K1", "K2", "KB", and "KW".
    If none, Ks are calculated with teh MyAMI model. Alternative Ks for non-
    seawater conditions are available in predefined NamedTuples. See file 
    "Constants" for details.
* pdict: dict
    Optional: can be used to provide some or all paramters as a dictionary/
    NamedTuples with the same key names. Any paramters in pdict will 
    overwrite manually specified paramters.

Returns
-------
NamedTuple containing all calculated parameters
"""
function boron_system(;
    pHtot=nothing, BT=nothing, BOH₃=nothing, BOH₄=nothing, ABT=nothing,
    ABOH₃=nothing, ABOH₄=nothing, δBT=nothing, δBOH₃=nothing, δBOH₄=nothing,
    alphaB=nothing, T_in=25.0, S_in=35.0, P_in =0.0, Ca=0.0102821, Mg=0.0528171, 
    ST=nothing, FT=nothing, pHsws=nothing, pHfree=nothing, 
    pHNBS=nothing, Ks=nothing, pdict=nothing, MyAMI_mode="calculate", 
    K_method="default", K_mode="static", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default", Ca_method="default", kwargs...)

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
        ST   = clean(isnothing(ST) ? calc_ST(S_in) : ST), # If not provided, calculated from S
        FT   = clean(isnothing(FT) ? calc_FT(S_in) : FT), # If not provided, calculated from S
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT, 
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        Mg = Mg,
        Ca = Ca,
        T_in = T_in,
        P_in = P_in,
        S_in = S_in,
        pHtot = pHtot,
        pHfree = pHfree,
        pHsws = pHsws,
        pHNBS = pHNBS
    )

    # If not provided, equilibrium constants are claculated with K calculator:
    if isnothing(Ks)

        new_Ks = K_calculator(; T_in, S_in, P_in=P_in, ST=nothing, FT=nothing,
        BT=nothing, K_method=K_method, KSO4_method=KSO4_method,
        BT_method=BT_method, KF_method=KF_method, KNH3_method=KNH3_method, 
        Ca_method=Ca_method, K_mode=K_mode, kwargs...)

        ps = merge(ps, new_Ks)

    else
        Ks = Ks isa Dict ? NamedTuple(Ks) : Ks
        ps = merge(ps, (; Ks=Ks))

    end


    # Calculate pH for all scales given an input pH value & scale:
    pH_results = calc_pH_scale(
        ps.pHtot, ps.pHfree, ps.pHsws, ps.pHNBS, ps.ST, ps.FT, ps.T_in .+ 273.15,
        ps.S_in, ps.Ks
    )
    if !isempty(pH_results)
        ps = merge(ps, pH_results)
    end

   # If pH is unknown, assign δBT value to calculate ABT, ABOH₃, and ABOH₄
    if isnothing(ps.pHtot)
        δBT = get(ps, :δBT, 39.61) # NEEDS CITATION 
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
        ps.pHtot, ps.pHfree, ps.pHsws, ps.pHNBS, ps.ST, ps.FT, ps.T_in .+ 273.15,
        ps.S_in, ps.Ks
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

    rem = (:pdict,) 
    ps = (; (k => v for (k, v) in pairs(ps) if k ∉ rem)...)

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
    to calculate MyAMI constants.
* T_in, S_in: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* P_in: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 
* Ks: NamedTuple
    Conatins named tuples of constants. Must contain:
    "K1", "K2", "KB", and "KW".
    If none, Ks are calculated with teh MyAMI model. Alternative Ks for non-
    seawater conditions are available in predefined NamedTuples. See file 
    "Constants" for details.
* pdict: dict
    Optional: can be used to provide some or all paramters as a dictionary/
    NamedTuples with the same key names. Any paramters in pdict will 
    overwrite manually specified paramters.

Returns
-------
NamedTuple containing all calculated parameters
"""
function boron_isotopes(;
    pHtot=nothing, BT=nothing, BOH₃=nothing, BOH₄=nothing, ABT=nothing, 
    ABOH₃=nothing, ABOH₄=nothing, δBT=nothing, δBOH₃=nothing, δBOH₄=nothing, 
    alphaB=nothing, T_in=25.0,  S_in=35.0, P_in=0.0, Ca=0.0102821, Mg=0.0528171,
    ST=nothing, FT=nothing, pHsws=nothing, pHfree=nothing, pHNBS=nothing, 
    Ks=nothing, pdict=nothing, MyAMI_mode="calculate", 
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default", Ca_method="default",
    K_mode="static", kwargs...
)

clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

ps = (
        ST   = clean(isnothing(ST) ? calc_ST(S_in) : ST), # If not provided, calculated from S
        FT   = clean(isnothing(FT) ? calc_FT(S_in) : FT), # If not provided, calculated from S
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
        T_in = T_in,
        P_in = P_in,
        S_in = S_in,
        pHtot = pHtot,
        pHfree = pHfree,
        pHsws = pHsws,
        pHNBS = pHNBS,
        kwargs...
    )

# If not provided, equilibrium constants are claculated with K calculator:
    if isnothing(Ks)

        new_Ks = K_calculator(; T_in, S_in, P_in=P_in, ST=nothing, FT=nothing,
        BT=nothing, K_method=K_method, KSO4_method=KSO4_method,
        BT_method=BT_method, KF_method=KF_method, KNH3_method=KNH3_method, 
        Ca_method=Ca_method, K_mode=K_mode, kwargs...)

        ps = merge(ps, new_Ks)

    else
        Ks = Ks isa Dict ? NamedTuple(Ks) : Ks
        ps = merge(ps, (; Ks=Ks))

    end

    # Calculate pH for all scales given an input pH value & scale:
    pH_results = calc_pH_scale(
        ps.pHtot, ps.pHfree, ps.pHsws, ps.pHNBS, ps.ST, ps.FT, ps.T_in + 273.15,
        ps.S_in, ps.Ks
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

    rem = (:pdict,) 
    ps = (; (k => v for (k, v) in pairs(ps) if k ∉ rem)...)

    return ps

end 


# Carbon & boron calculations
"""
Calculates the carbon and boron species as well as boron isotopes of seawater 
from given parameters. 
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
    to calculate MyAMI constants.
* T_in, S_in: array-like
    Temperature in Celcius and salinity in PSU for the condtions that the
    measurments were taken in. Used in calculating MyAMI constants.
* P_in: array-like
    Pressure in Bar for the conditions that the measuremnts were taken in. Used
    in pressure-correcting constants. 
* T_out, S_out: array-like
    Temperature in Celcius and salinity in PSU of desired output conditions. Used
    in calculating constants.
* P_out: array-like
    Pressure in bar of desired output conditions. Used in pressure-correcting
    constants.
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
* pdict: dict
    Optional: can be used to provide some or all paramters as a dictionary/
    NamedTuples with the same key names. Any paramters in pdict will 
    overwrite manually specified paramters.

Returns
-------
NamedTuple containing all calculated parameters

"""
function whole_system(;
    pH=nothing, pHtot=nothing, DIC=nothing, TA=nothing, CO₂=nothing, HCO₃=nothing,
    CO₃=nothing, pCO₂=nothing, fCO₂=nothing, BT =nothing, BOH₃=nothing,
    BOH₄=nothing, ABT=nothing, ABOH₃=nothing, ABOH₄=nothing, δBT=nothing,
    δBOH₃=nothing, δBOH₄=nothing, alphaB=nothing, Ca=0.0102821,
    Mg=0.0528171, T_in=25.0, T_out=nothing, S_in=35.0, S_out=nothing,
    P_in=0.0, P_out=nothing, PT=0.0, SiT=0.0, H2ST=0.0, NH4T=0.0, ST=nothing, FT=nothing,
    pHsws=nothing, pHfree=nothing, pHNBS=nothing, unit="umol", scale="total", Ks=nothing,
    pdict=nothing, ΩC=nothing, ΩA=nothing, MyAMI_mode="calculate",
    K_method="default", KSO4_method="default", BT_method="default", 
    KF_method="default", KNH3_method="default", Ca_method="default",
    K_mode="static", kwargs...
)

    udict = Dict(
        "mol" => 1.0,
        "mmol" => 1.0e3,
        "umol" => 1.0e6,
        "nmol" => 1.0e9,
        "pmol" => 1.0e12,
        "fmol" => 1.0e15
    )

    m = get(udict, unit, 1.0)
    scale_unit(v) = isnothing(v) ? nothing : v / m
    scale_gas(v) = isnothing(v) ? nothing : v / 1e6 
    clean(x) = isnothing(x) ? nothing : (x < 0.0 ? NaN : x)

    if isnothing(Ks)
        Ks_env = K_calculator(;
            T_in=T_in, S_in=S_in, P_in=P_in, ST=ST, FT=FT, BT=BT,
            K_method=K_method, KSO4_method=KSO4_method, BT_method=BT_method,
            KF_method=KF_method, KNH3_method=KNH3_method, Ca_method=Ca_method, K_mode=K_mode
        )
    else
        Ks_env = Ks
    end

    ST_val = Ks_env.ST
    FT_val = Ks_env.FT
    KS_val = Ks_env.Ks.KS
    KF_val = Ks_env.Ks.KF
    
    BT_from_env = hasproperty(Ks_env, :BT) ? Ks_env.BT : 
                  (hasproperty(Ks_env.Ks, :BT) ? Ks_env.Ks.BT : 0.0)

    internal_BT = !isnothing(BT) ? clean(scale_unit(BT)) : 
                  (BT_from_env > 1.0 ? BT_from_env / 1e6 : BT_from_env)

    fH_val = Ks_env.fH isa AbstractArray ? Ks_env.fH[1] : Ks_env.fH
    
    tf_fac = (1.0 + ST_val / KS_val) 
    ts_fac = tf_fac / (1.0 + ST_val / KS_val + FT_val / KF_val)
    
    # --- THE WAY IN: Map generic pH, then strictly convert to pHtot ---
    if !isnothing(pH) && isnothing(pHtot) && isnothing(pHsws) && isnothing(pHfree) && isnothing(pHNBS)
        if lowercase(scale) == "total"
            pHtot = pH
        elseif lowercase(scale) == "sws"
            pHsws = pH
        elseif lowercase(scale) == "free"
            pHfree = pH
        elseif lowercase(scale) == "nbs"
            pHNBS = pH
        end
    end

    if isnothing(pHtot)
        if !isnothing(pHsws)
            pHtot = pHsws - log10(ts_fac)
        elseif !isnothing(pHfree)
            pHtot = pHfree - log10(tf_fac)
        elseif !isnothing(pHNBS)
            pHtot = pHNBS + log10(fH_val) - log10(ts_fac)
        end
    end
    
    local_H = !isnothing(pHtot) ? 10.0^(-pHtot) : nothing

    ps = (;
        kwargs...,
        DIC = clean(scale_unit(DIC)),
        TA = scale_unit(TA),
        CO₂ = clean(scale_unit(CO₂)),
        HCO₃ = clean(scale_unit(HCO₃)),
        CO₃ = clean(scale_unit(CO₃)),
        PT = clean(scale_unit(PT)),
        SiT = clean(scale_unit(SiT)),
        pCO₂ = clean(scale_gas(pCO₂)),
        fCO₂ = clean(scale_gas(fCO₂)),
        BT = internal_BT,
        BOH₃ = clean(scale_unit(BOH₃)),
        BOH₄ = clean(scale_unit(BOH₄)),
        ST = Ks_env.ST,
        FT = Ks_env.FT,
        Mg = Mg,
        Ca = Ca,
        fH = fH_val,
        δBT = δBT,
        δBOH₃ = δBOH₃,
        δBOH₄ = δBOH₄,
        ABT = ABT, 
        ABOH₃ = ABOH₃,
        ABOH₄ = ABOH₄,
        alphaB = alphaB,
        H2ST = clean(scale_unit(H2ST)),
        NH4T = clean(scale_unit(NH4T)),
        T_in = T_in,
        P_in = P_in,
        S_in = S_in,
        T_out = T_out,
        S_out = S_out,
        P_out = P_out,
        pH = nothing,     # <-- Added to explicitly clear generic pH
        pHtot = pHtot,
        pHfree = nothing, # Cleared to prevent inner solvers from touching them
        pHsws = nothing,
        pHNBS = nothing,
        scale = "total", # Force the internal scale to be total
        unit = unit,
        Ks = Ks_env.Ks
    )

    alphaB   = !isnothing(get(ps, :alphaB, nothing))   ? ps.alphaB   : Isotopes.get_alphaB()
    δBT_val = (isnothing(δBT) && isnothing(ABT)) ? 39.61 : δBT
    ABT_val   = !isnothing(δBT_val)   ? Isotopes.δ11_to_A11(δBT_val)   : ABT
    ABOH₃_val = !isnothing(δBOH₃) ? Isotopes.δ11_to_A11(δBOH₃) : ABOH₃
    ABOH₄_val = !isnothing(δBOH₄) ? Isotopes.δ11_to_A11(δBOH₄) : ABOH₄

    ps = merge(ps, (δBT = δBT_val, ABT = ABT_val, ABOH₃ = ABOH₃_val, ABOH₄ = ABOH₄_val))
    nBiso = count(!isnothing, (ABT,)) + count(!isnothing, (ABOH₃, ABOH₄))

    local_fCO2 = ps.fCO₂
    local_pCO2 = ps.pCO₂
    local_CO2  = ps.CO₂

    if isnothing(local_CO2)
        if !isnothing(local_fCO2)
            local_CO2  = local_fCO2 * ps.Ks.K0
            local_pCO2 = isnothing(local_pCO2) ? fCO₂_to_pCO₂(local_fCO2, ps.T_in) : local_pCO2
        elseif !isnothing(local_pCO2)
            local_fCO2 = pCO₂_to_fCO₂(local_pCO2, ps.T_in)
            local_CO2  = local_fCO2 * ps.Ks.K0
            local_pCO2 = local_pCO2
        end
    end

    ps = merge(ps, (CO₂ = local_CO2, fCO₂ = local_fCO2, pCO₂ = local_pCO2))

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

    # --- THE WAY OUT: Calculate other pH scales from the finished Total pH ---
    solved_pHtot = ps.pHtot
    final_pHsws  = solved_pHtot + log10(ts_fac)
    final_pHfree = solved_pHtot + log10(tf_fac)
    final_pHNBS  = final_pHsws - log10(fH_val) 

    ps = merge(ps, (;
        pHfree = final_pHfree,
        pHsws  = final_pHsws,
        pHNBS  = final_pHNBS
    ))

    final_δBT   = !isnothing(δBT) ? δBT : A11_to_δ11(ps.ABT)
    final_δBOH₃ = !isnothing(δBOH₃) ? δBOH₃ : A11_to_δ11(ps.ABOH₃)
    final_δBOH₄ = !isnothing(δBOH₄) ? δBOH₄ : A11_to_δ11(ps.ABOH₄)

    ps = merge(ps, (; δBT = final_δBT, δBOH₃ = final_δBOH₃, δBOH₄ = final_δBOH₄))

    rf = calc_revelle_factor(ps.TA, ps.DIC, ps.BT, ps.PT, ps.SiT, ps.ST, 
                             ps.FT, ps.H2ST, ps.NH4T, ps.Ks)
    ps = merge(ps, (revelle_factor = rf,))

    oCa = if !isnothing(get(ps, :Ca, nothing))
        ps.Ca * ps.S_in / 35.0
    else
        0.0102821 * ps.S_in / 35.0
    end

    ΩA = ps.CO₃ * oCa / ps.Ks.KspA
    ΩC = ps.CO₃ * oCa / ps.Ks.KspC
    ps = merge(ps, (ΩA = ΩA, ΩC = ΩC))

    if m ≠ 1
        target_keys = (:DIC, :TA, :BT, :CO₂, :HCO₃, :CO₃, :PT, :SiT, :BOH₃, :BOH₄, :CAlk, :BAlk, :PAlk, :OH, :SiAlk, :HSO₄, :Hfree, :HF)
        
        rescaled_concs = (; [
            begin
                val = getfield(ps, k)
                if isnothing(val)
                    k => nothing
                elseif (k === :BT || k === :DIC || k === :TA) && val > 1.0
                    k => val
                else
                    k => val * m
                end
            end for k in target_keys if hasproperty(ps, k)
        ]...)
        
        ps = merge(ps, rescaled_concs)
    end

    rescaled_gases = (
        pCO₂ = ps.pCO₂ * 1e6,
        fCO₂ = ps.fCO₂ * 1e6
    )
    ps = merge(ps, rescaled_gases)

    if (!isnothing(ps.T_out)) || (!isnothing(ps.S_out)) || (!isnothing(ps.P_out))
        out_params = (
            T_out = isnothing(ps.T_out) ? ps.T_in : ps.T_out,
            S_out = isnothing(ps.S_out) ? ps.S_in : ps.S_out,
            P_out = isnothing(ps.P_out) ? ps.P_in : ps.P_out
        )
        ps = merge(ps, out_params)

        local_BT = ps.S_in ≠ ps.S_out ? ps.BT .* ps.S_out ./ ps.S_in : ps.BT
        local_ST = ps.S_in ≠ ps.S_out ? ps.ST .* ps.S_out ./ ps.S_in : ps.ST
        local_FT = ps.S_in ≠ ps.S_out ? ps.FT .* ps.S_out ./ ps.S_in : ps.FT
    
        out_cond = whole_system(;
            TA = ps.TA,
            DIC = ps.DIC,
            T_in = ps.T_out,
            S_in = ps.S_out,
            P_in = ps.P_out,
            T_out = nothing, 
            S_out = nothing,
            P_out = nothing,
            unit = ps.unit, #"mol",
            scale = ps.scale,
            Ca = ps.Ca,
            Mg = ps.Mg,
            BT = local_BT,
            FT = local_FT,
            ST = local_ST,
            PT = ps.PT,
            SiT = ps.SiT,
            MyAMI_mode = MyAMI_mode,
            K_method="default",
            KSO4_method="default",
            BT_method="default"
        )

        outputs = [
            "BAlk", "BT", "CAlk", "CO₂", "CO₃", "DIC", "H", "HCO₃", "HF", "HSO₄", 
            "Hfree", "Ks", "OH", "PAlk", "SiAlk", "TA", "FT","PT", "ST", "SiT", 
            "fCO₂", "pCO₂", "pHfree", "pHsws", "pHtot", "pHNBS", "ΩA", "ΩC", 
            "revelle_factor", "BOH₃", "BOH₄", "δBT", "δBOH₃", "δBOH₄", "ABT",
            "ABOH₃", "ABOH₄", "alphaB"
        ]

        in_vals = (; (Symbol(k * "_in") => getfield(ps, Symbol(k)) 
                      for k in outputs if hasproperty(ps, Symbol(k)))...)

        out_vals = (; (Symbol(k) => getfield(out_cond, Symbol(k)) 
                       for k in outputs if hasproperty(out_cond, Symbol(k)))...)
        ps = merge(ps, in_vals, out_vals)
    end 
    
    keys_to_keep = Tuple(k for k in keys(ps) if k ∉ (:pdict, :unit, :scale))
    ps = NamedTuple{keys_to_keep}(ps)
    
    return ps
end

end # module