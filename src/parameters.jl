# What the inputs are: every parameter the package accepts, and how they relate.
#
# One file, because these tables are read by three things that must agree — the solver's
# accepted names and defaults, the input validation, and the solve-order dispatch. Split
# across files, their descriptions of "what determines the system" drift apart.

"""
Every parameter a calculation accepts, with its default.

The single declaration of the input surface: the solver's accepted names, its defaults, and
what the validation checks against all read from here.
"""
const PARAMETER_DEFAULTS = (
    # Carbonate system parameters.
    TA = nothing, DIC = nothing,
    pHtot = nothing, pHsws = nothing, pHfree = nothing, pHNBS = nothing,
    CO₂ = nothing, HCO₃ = nothing, CO₃ = nothing, pCO₂ = nothing, fCO₂ = nothing,
    ΩA = nothing, ΩC = nothing,

    # Totals and seawater composition.
    BT = nothing, ST = nothing, FT = nothing,
    PT = 0.0, SiT = 0.0, H2ST = 0.0, NH4T = 0.0,
    Ca = nothing, Mg = nothing,

    # Conditions.
    temp_c = 25.0, sal = 35.0, pres_bar = 0.0,

    # How the constants are calculated, and what units the values are in.
    unit = "umol", Ks = nothing, MyAMI_mode = "approximate",
    K_method = "default", KSO4_method = "default", BT_method = "default",
    KF_method = "default", KNH3_method = "default", Ca_method = "default",

    # Boron speciation.
    BOH₃ = nothing, BOH₄ = nothing,

    # Boron isotopes.
    δBT = nothing, δBOH₃ = nothing, δBOH₄ = nothing,
    ABT = nothing, ABOH₃ = nothing, ABOH₄ = nothing, alphaB = nothing,
)

"""
The scope each parameter belongs to. Anything absent applies to every scope.

A solver only accepts the parameters its scope covers, so asking `carbon_system` for `δBT`
is an unrecognised argument rather than a value silently ignored — and equally, a boron-only
solver does not accept `DIC`.

Absent, and so shared by every scope: pH on any scale (the unknown all three subsystems
share), the totals `BT`/`ST`/`FT`, the seawater composition `Ca`/`Mg`, the conditions, and
every setting. Those describe the water rather than one subsystem's chemistry.
"""
const PARAMETER_SCOPE = (
    # Carbonate system, including the nutrient contributions to alkalinity.
    TA = :carbon, DIC = :carbon, CO₂ = :carbon, HCO₃ = :carbon, CO₃ = :carbon,
    pCO₂ = :carbon, fCO₂ = :carbon, ΩA = :carbon, ΩC = :carbon,
    PT = :carbon, SiT = :carbon, H2ST = :carbon, NH4T = :carbon,
    # Boron speciation.
    BOH₃ = :boron, BOH₄ = :boron,
    # Boron isotopes.
    δBT = :isotopes, δBOH₃ = :isotopes, δBOH₄ = :isotopes,
    ABT = :isotopes, ABOH₃ = :isotopes, ABOH₄ = :isotopes, alphaB = :isotopes,
)

"How many independent constraints a subsystem needs before it can be solved for pH."
const CONSTRAINTS_NEEDED = (carbon = 2, boron = 1, isotopes = 1)

"The parameter names a solver of this scope accepts."
function _accepted_parameters(scope::Tuple{Vararg{Symbol}})
    return Tuple(name for name in keys(PARAMETER_DEFAULTS)
                 if !haskey(PARAMETER_SCOPE, name) ||
                    getproperty(PARAMETER_SCOPE, name) in scope)
end

"""
Parameters grouped by the degree of freedom each constrains.

Two members of one group are not two measurements. pH on four scales is one measurement
expressed four ways; ΩA and ΩC are both statements about [CO₃²⁻]; pCO₂ and fCO₂ are both
statements about dissolved CO₂; δ and A are the same isotope value in two notations.

**`BT`, `δBT` and `ABT` are deliberately absent.** They are totals with defaults — `BT` from
salinity, `δBT` from modern seawater — so supplying one constrains nothing on its own. It is
the *speciated* member that carries the information.
"""
const PARAMETER_GROUPS = (
    pH    = (:pHtot, :pHsws, :pHfree, :pHNBS),
    CO₂   = (:CO₂, :pCO₂, :fCO₂),
    CO₃   = (:CO₃, :ΩA, :ΩC),
    HCO₃  = (:HCO₃,),
    DIC   = (:DIC,),
    TA    = (:TA,),
    BOH₃  = (:BOH₃,),
    BOH₄  = (:BOH₄,),
    δBOH₃ = (:δBOH₃, :ABOH₃),
    δBOH₄ = (:δBOH₄, :ABOH₄),
)

"""
Which subsystem each group constrains.

`pH` is `:shared` because it is the unknown that links the three: fix it in any one
subsystem and the other two follow. That is the whole basis of the solve order below.
"""
const GROUP_SUBSYSTEM = (pH = :shared, CO₂ = :carbon, CO₃ = :carbon, HCO₃ = :carbon,
                         DIC = :carbon, TA = :carbon, BOH₃ = :boron, BOH₄ = :boron,
                         δBOH₃ = :isotopes, δBOH₄ = :isotopes)

"""
Concentrations held internally in mol/kg that must be converted back to the caller's unit on
the way out, by scope.

`BT` is shared because every scope carries total boron — the carbonate system needs it for
alkalinity, and the boron system is about it.
"""
const SHARED_CONCENTRATIONS = (:BT,)

const CARBON_CONCENTRATIONS = (:DIC, :TA, :CO₂, :HCO₃, :CO₃, :PT, :SiT,
                               :CAlk, :BAlk, :PAlk, :OH, :SiAlk, :HSO₄, :Hfree, :HF,
                               :Alk_H2S, :Alk_NH3)

const BORON_CONCENTRATIONS = (:BOH₃, :BOH₄)

