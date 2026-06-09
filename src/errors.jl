# src/errors.jl
using ForwardDiff, LinearAlgebra

"""
    propagate_errors(target_func; inputs::NamedTuple, errors::NamedTuple)

Propagates uncertainties through a carbonate chemistry function using 
First-Order Taylor Series expansion (Automatic Differentiation).
"""
function propagate_errors(target_func; inputs::NamedTuple, errors::NamedTuple)
    error_keys = keys(errors)
    # Ensure we are working with Floats for the AD process
    base_values = Float64[getproperty(inputs, k) for k in error_keys]
    uncertainties = Float64[getproperty(errors, k) for k in error_keys]

    for key in keys(errors)
    if getproperty(inputs, key) === nothing
        error("Parameter '$key' was given an error/uncertainty, but its input value is 'nothing'. " *
              "Please provide a numerical value for '$key' in your inputs.")
    end
end

    # 1. Define the AD-compatible wrapper
    function wrapped_math(x_vec)
        # Create a mutable copy of the inputs
        perturbed_inputs = Dict{Symbol, Any}(pairs(inputs))
        
        # Inject the "Dual Numbers" from ForwardDiff
        for (i, key) in enumerate(error_keys)
            perturbed_inputs[key] = x_vec[i]
        end
        
        # Run the calculation
        res = target_func(; perturbed_inputs...)
        
        # ForwardDiff needs a Vector of numbers as output.
        # We filter for only the numeric results (ignoring things like the Ks struct)
        return [v for v in values(res) if v isa Number]
    end

    # 2. Calculate the Jacobian (The matrix of partial derivatives)
    jac = ForwardDiff.jacobian(wrapped_math, base_values)

    # 3. Perform the propagation: σ_out = sqrt( Σ (∂out/∂in * σ_in)^2 )
    out_variances = (jac .^ 2) * (uncertainties .^ 2)
    # Sum across the rows (the input contributions) and take sqrt
    out_errors_vec = sqrt.(sum(out_variances, dims=2))

    # 4. Reconstruct the results as a NamedTuple
    baseline_res = target_func(; inputs...)
    res_keys = [k for (k, v) in pairs(baseline_res) if v isa Number]
    
    err_tuple = NamedTuple{Tuple(res_keys)}(vec(out_errors_vec))

    return (val = baseline_res, err = err_tuple)
end