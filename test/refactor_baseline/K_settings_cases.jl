# The cases behind the K_settings baseline, shared by the script that captures it and the test
# that checks against it so the two cannot drift apart.
#
# Written in the *old* API only (flat `K_method = …` keywords) on purpose: the capture script
# has to run on the code from before `K_settings` existed, and the check then runs the same
# calls on the new code — once as written, and once translated into `K_settings` by hand.
#
# Each case is `name => (kind, kwargs)`:
#   :solve      carbon_system(; kwargs...), recording every number in `.val`
#   :constants  calculate_constants(; kwargs...), recording every constant and total

const OLD_K_METHODS = ["default", "KGen", "Roy 1993", "GP 1989", "Hansson 1973", "DM 1987",
                       "HM 1973", "Mehrbach 1973 A", "Mehrbach 1973 B", "CW 2003",
                       "Lueker 2000", "MPM 2002", "Millero 2002", "Millero 2006",
                       "Millero 2010", "Waters 2014", "SB 2020", "Papadimitriou 2018",
                       "Sulpis 2020"]

function K_settings_cases()
    cases = Pair{String, Tuple{Symbol, NamedTuple}}[]
    add!(name, kind, kw) = push!(cases, name => (kind, kw))
    pair = (TA = 2300.0, DIC = 2000.0)

    # Every K1/K2 parameterisation over a small grid of conditions. Millero 1979 is pure
    # water only, so it gets its own grid at salinity 0.
    for k in OLD_K_METHODS, temp_c in (2.0, 25.0), sal in (25.0, 35.0), pres_bar in (0.0, 300.0)
        conditions = (; temp_c, sal, pres_bar, K_method = k)
        tag = "$(k)_T$(temp_c)_S$(sal)_P$(pres_bar)"
        add!("solve_$tag", :solve, merge(pair, conditions))
        add!("constants_$tag", :constants, conditions)
    end
    for temp_c in (2.0, 25.0), pres_bar in (0.0, 300.0)
        conditions = (; temp_c, sal = 0.0, pres_bar, K_method = "Millero 1979")
        tag = "Millero 1979_T$(temp_c)_S0.0_P$(pres_bar)"
        add!("solve_$tag", :solve, merge(pair, conditions))
        add!("constants_$tag", :constants, conditions)
    end

    # pH on the NBS scale goes through fH, which Mehrbach 1973 B and Millero 1979 replaced.
    for k in ("default", "Lueker 2000", "Mehrbach 1973 A", "Mehrbach 1973 B")
        add!("nbs_$k", :solve, (TA = 2300.0, pHNBS = 8.1, temp_c = 15.0, K_method = k))
    end
    add!("nbs_Millero 1979", :solve,
         (TA = 2300.0, pHNBS = 8.1, temp_c = 15.0, sal = 0.0, K_method = "Millero 1979"))

    # The other choices, each crossed with the ones it could interact with.
    for kso4 in ("default", "Dickson", "Khoo", "WM13"), kf in ("default", "Dickson", "Perez"),
        bt in ("default", "Uppstrom", "Lee", "KSK18"), pres_bar in (0.0, 300.0)
        kw = (; temp_c = 15.0, sal = 35.0, pres_bar, K_method = "Lueker 2000", KSO4_method = kso4,
              KF_method = kf, BT_method = bt)
        tag = "$(kso4)_$(kf)_$(bt)_P$(pres_bar)"
        add!("solve_Lueker_$tag", :solve, merge(pair, kw))
        add!("constants_Lueker_$tag", :constants, kw)
    end
    for knh3 in ("default", "Millero", "Clegg"), ca in ("default", "modern", "Culkin", "RT67"),
        k in ("default", "Lueker 2000")
        kw = (; temp_c = 15.0, sal = 33.0, K_method = k, KNH3_method = knh3, Ca_method = ca)
        tag = "$(k)_$(knh3)_$(ca)"
        add!("solve_nh3_$tag", :solve, merge(pair, kw, (NH4T = 5.0,)))
        add!("constants_nh3_$tag", :constants, kw)
    end

    # Nutrients reach PAlk and SiAlk through constants Millero 1979 used to zero, so it is
    # deliberately absent: with PT > 0 the old zeroed constants gave PAlk = -PT.
    for k in OLD_K_METHODS
        add!("nutrients_$k", :solve, merge(pair, (temp_c = 15.0, sal = 35.0, PT = 2.0,
                                                  SiT = 20.0, H2ST = 1.0, NH4T = 3.0,
                                                  K_method = k)))
    end

    # Overrides that used to win: Mehrbach 1973 A's Culkin boron over BT_method, and a
    # supplied BT over both.
    for k in ("Mehrbach 1973 A", "Mehrbach 1973 B"), bt in ("Lee", "KSK18")
        add!("override_$(k)_$bt", :solve, merge(pair, (temp_c = 15.0, K_method = k,
                                                       BT_method = bt)))
    end
    add!("override_supplied_BT", :solve, merge(pair, (temp_c = 15.0, BT = 400.0,
                                                      K_method = "Mehrbach 1973 A")))

    # Composition: Kgen's correction, and an explicit Ca with a fitted parameterisation.
    add!("composition_KGen", :solve, merge(pair, (Ca = 0.012, Mg = 0.05)))
    add!("composition_KGen_approximate", :solve,
         merge(pair, (Ca = 0.012, Mg = 0.05, K_method = "KGen", MyAMI_mode = "approximate")))
    add!("composition_Lueker", :solve,
         merge(pair, (Ca = 0.012, K_method = "Lueker 2000")))

    return cases
end

"Run one case, returning `field => value` for every number it produces."
function run_K_settings_case(kind::Symbol, kw::NamedTuple)
    if kind === :solve
        val = carbon_system(; kw...).val
        return [String(k) => Float64(v) for (k, v) in pairs(val) if v isa Real]
    else
        env = CarbonateCalculator.calculate_constants(; kw...)
        fields = [String(k) => Float64(v) for (k, v) in pairs(env.Ks)]
        for k in (:BT, :ST, :FT, :Ca, :Mg, :fH, :tf_fac, :ts_fac)
            push!(fields, String(k) => Float64(getproperty(env, k)))
        end
        return fields
    end
end
