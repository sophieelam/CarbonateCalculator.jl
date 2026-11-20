# Calcium and saturation state, checked against the published formulas directly.
#
# Saturation state is CO₃²⁻ × Ca²⁺ / Ksp. The CO2SYS comparisons cannot isolate the Ca²⁺
# part: CO2SYS uses Culkin calcium only with its GEOSECS options, which differ from this
# package in several other documented ways, so a calcium error there hides among them. That
# is how Culkin and RT67 came to be diluted twice (Ω low by S/35) without a failing test.
# Here each option's in-situ calcium is written out from its reference, independently of the
# package, and Ω must follow from it at salinities away from 35.

using Test, CarbonateCalculator

"In-situ calcium, mol/kg-SW, as each reference writes it."
const CALCIUM_FORMULAS = Dict(
    "modern" => S -> CarbonateCalculator.MODERN_CALCIUM * S / 35,
    "Culkin" => S -> 0.01026 * S / 35,                    # Culkin (1965); CO2SYS options 6, 7
    "RT67"   => S -> 0.02128 / 40.087 * S / 1.80655,      # Riley & Tongudai (1967); CO2SYS
    "CO2SYS" => S -> 0.02128 / 40.087 * S / 1.80655,
)

@testset "Calcium and saturation state" begin
    for (method, formula) in CALCIUM_FORMULAS, sal in (5.0, 20.0, 35.0, 42.0)
        r = carbon_system(; TA = 2300.0, DIC = 2000.0, temp_c = 15.0, sal,
                            K_settings = (K1_K2_method = "Lueker 2000", Ca_method = method))
        CO3 = r.val.CO₃ * 1e-6                       # µmol/kg → mol/kg
        Ca = formula(sal)
        @test r.val.ΩC ≈ CO3 * Ca / r.Ks.KspC rtol = 1e-12
        @test r.val.ΩA ≈ CO3 * Ca / r.Ks.KspA rtol = 1e-12
        # Reported calcium is at S = 35, as a supplied `Ca` is (README, Calcium Concentration).
        @test r.val.Ca ≈ formula(35.0) rtol = 1e-12
    end

    # "CO2SYS" is RT67 under another name: identical, not merely close.
    a = carbon_system(; TA = 2300.0, DIC = 2000.0, temp_c = 15.0, sal = 30.0,
                        K_settings = (Ca_method = "CO2SYS",))
    b = carbon_system(; TA = 2300.0, DIC = 2000.0, temp_c = 15.0, sal = 30.0,
                        K_settings = (Ca_method = "RT67",))
    @test a.val.ΩA == b.val.ΩA

    # A supplied Ca is at S = 35 and wins over Ca_method.
    r = carbon_system(; TA = 2300.0, DIC = 2000.0, temp_c = 15.0, sal = 20.0, Ca = 0.0103,
                        K_settings = (K1_K2_method = "Lueker 2000", Ca_method = "Culkin"))
    @test r.val.ΩC ≈ r.val.CO₃ * 1e-6 * (0.0103 * 20 / 35) / r.Ks.KspC rtol = 1e-12
end
