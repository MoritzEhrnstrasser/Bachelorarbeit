using CSV, DataFrames
using Dates
using Statistics
using LinearAlgebra
using Distributions

function matrix_steps_timed(X, y, p)
    N = length(y)
    Xmat = hcat(ones(size(X, 1)), X)
    Y = y

    t1 = time()
    XtX_inv = inv(transpose(Xmat) * Xmat)
    beta_hat = XtX_inv * transpose(Xmat) * y
    t2 = time()

    y_hat = Xmat * beta_hat
    residuals = Y .- y_hat
    t3 = time()

    sigma2_hat = sum(residuals .^ 2) / (N - 2)
    se_beta = sqrt.(sigma2_hat .* diag(XtX_inv))
    t4 = time()

    t_stats = beta_hat ./ se_beta
    p_values = 2 .* (1 .- cdf(TDist(N - p - 1), abs.(t_stats)))
    t5 = time()

    return (
        inversion = t2 - t1,
        prediction = t3 - t2,
        standard_error = t4 - t3,
        p_values = t5 - t4
    )
end

data_path = "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"
files = filter(f -> endswith(f, ".csv"), readdir(data_path; join=true))

results = DataFrame(N=Int[], p=Int[],
                    inversion=Float64[], prediction=Float64[],
                    standard_error=Float64[], p_values=Float64[])

for file in files
    df = CSV.read(file, DataFrame)
    y = df.y
    X = Matrix(df[:, occursin.("x", names(df))])

    m = match(r"N(\d+)_p(\d+)", file)
    if m === nothing
        @warn "Dateiname entspricht nicht dem erwarteten Format: $file"
        continue  # Datei überspringen
    end
    N = parse(Int, m.captures[1])
    p = parse(Int, m.captures[2])
    

    step_results = [matrix_steps_timed(X, y, p) for _ in 1:500]
    tmp_df = DataFrame(step_results)
    means = [mean(tmp_df[!, col]) for col in names(tmp_df)]
    means = means * 1e6 # Convert to microseconds
    means = round.(means, digits=2)

    push!(results, (N, p, means...))
end

CSV.write(joinpath(data_path, "timing_steps_julia.csv"), results)
