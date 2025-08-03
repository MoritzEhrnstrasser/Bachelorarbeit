

using CSV, DataFrames, LinearAlgebra, Statistics, Glob

# 1. Konfiguration
data_path = "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"
files = glob("*.csv", data_path)
isample = 500

# 2. Regressionsfunktion mit interner Zeitmessung (via time())
function matrix_regression_timed(X::Matrix{Float64}, y::Vector{Float64})
    
    N = size(X, 1)
    p = size(X, 2)
    Xmat = hcat(ones(N), X)

    start_time = time()

    XtX_inv = inv(Xmat' * Xmat)
    β = XtX_inv * (Xmat' * y)

    ŷ = Xmat * β
    residuals = y - ŷ
    σ² = sum(residuals .^ 2) / (N - p - 1)
    se = sqrt.(diag(σ² * XtX_inv))

    t_stats = β ./ se
    p_values = 2 .* (1 .- cdf(TDist(N - p - 1), abs.(t_stats)))

    end_time = time()
    return end_time - start_time  # Dauer in Sekunden
end

# 3. DataFrame für Ergebnisse
results = DataFrame(N = Int[], p = Int[], median_time_us = Float64[], sd_time_us = Float64[])

# 4. Loop über alle Dateien
for file in files
    # Metadaten aus Dateiname
    caps = match(r"N(\d+)_p(\d+)", basename(file)).captures
    N = parse(Int, caps[1])
    p = parse(Int, caps[2])

    # Daten einlesen
    df = CSV.read(file, DataFrame)
    y = df.y
    X = Matrix(df[:, r"^x"])

    # Mehrfache Zeitmessung
    times = [ matrix_regression_timed(X, y) for _ in 1:isample ]

    # Umrechnen in Mikrosekunden
    median_us = median(times) * 1e6
    sd_us     = std(times)    * 1e6

    push!(results, (N=N, p=p, median_time_us=median_us, sd_time_us=sd_us))
end

# 5. Speichern
CSV.write(joinpath(data_path, "julia_whole_timed.csv"), results)

