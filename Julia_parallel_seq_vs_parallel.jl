# ==== Julia-Skript: Seq. vs. Parallel mit Distributed ====

using CSV, DataFrames, LinearAlgebra, Statistics, Distributions, Distributed

# 1. Setup
data_path = "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"
files = readdir(data_path; join=true) |> x -> filter(f->endswith(f, ".csv"), x)
isample = 3000
n_cores = 4

# 2. Starte Worker
addprocs(n_cores)
@everywhere using CSV, DataFrames, LinearAlgebra, Statistics, Distributions

# 3. Regressions- & Benchmark-Funktion
@everywhere function matrix_regression_timed(X::Matrix{Float64}, y::Vector{Float64})
    t1 = time()
    N, p = size(X)
    Xmat = hcat(ones(N), X)
    XtX_inv = inv(Xmat' * Xmat)
    β = XtX_inv * (Xmat' * y)
    ŷ = Xmat * β
    residuals = y - ŷ
    σ² = sum(residuals .^ 2) / (N - p - 1)
    se = sqrt.(diag(σ² * XtX_inv))
    t_stats = β ./ se
    p_values = 2 .* (1 .- cdf(TDist(N - p - 1), abs.(t_stats)))
    return time() - t1
end

@everywhere function benchmark_single_file(file::String, isample::Int)
    df = CSV.read(file, DataFrame)
    y  = df.y
    X  = Matrix(df[:, r"^x"])
    caps = match(r"N(\d+)_p(\d+)", basename(file)).captures
    N    = parse(Int, caps[1])
    p    = parse(Int, caps[2])
    times = [matrix_regression_timed(X, y) for _=1:isample]
    DataFrame(
      N = N,
      p = p,
      median_time_us = median(times) * 1e6,
      sd_time_us     = std(times)    * 1e6
    )
end

# 4.a Sequentiell
t_start_seq = time()
results_seq = [benchmark_single_file(f, isample) for f in files]
t_end_seq = time()
t_seq = t_end_seq - t_start_seq

# 4.b Parallel mit pmap
t_start_par = time()
results_par = pmap(f->benchmark_single_file(f, isample), files)
t_end_par = time()
t_par = t_end_par - t_start_par

# 5. Zusammenführen
results_seq_df = vcat(results_seq...)
results_par_df = vcat(results_par...)

# 6. Speedup & Effizienz
speedup    = t_seq / t_par
efficiency = speedup / n_cores

println("Julia Benchmark (n_cores = $n_cores)")
println(" Sequentiell: $(round(t_seq, digits=2)) s")
println(" Parallel:    $(round(t_par, digits=2)) s")
println(" Speedup:     $(round(speedup, digits=2))")
println(" Effizienz:   $(round(efficiency*100, digits=1)) %")

# 7. Export optional
CSV.write("/tmp/timing_seq_julia.csv", results_seq_df)
CSV.write("/tmp/timing_par_julia.csv", results_par_df)
