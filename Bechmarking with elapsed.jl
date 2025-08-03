using CSV, DataFrames, LinearAlgebra, Statistics, Glob

# 1. Konfiguration
data_path = "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"
files = glob("*.csv", data_path)
isample = 500

# 2. Matrixbasierte Regressionsfunktion (ohne Zeitmessung)
function matrix_regression!(X::Matrix{Float64}, y::Vector{Float64})
    N = size(X, 1)
    Xmat = hcat(ones(N), X)
    β = (Xmat' * Xmat) \ (Xmat' * y)
    ŷ = Xmat * β
    residuals = y - ŷ
    σ² = sum(residuals .^ 2) / (N - 2)
    XtX_inv = inv(Xmat' * Xmat)
    se = sqrt.(diag(σ² * XtX_inv))
    # t_stats und p_values berechnen, aber nicht zurückgeben
    t_stats = β ./ se
    p_values = 2 .* (1 .- cdf(TDist(N - p - 1), abs.(t_stats)))
    return nothing
end

# 3. DataFrame für Ergebnisse
results = DataFrame(N = Int[], p = Int[], median_time_us = Float64[], sd_time_us = Float64[])

# 4. Loop über alle Dateien
for file in files
    # Metadaten aus dem Dateinamen extrahieren
    caps = match(r"N(\d+)_p(\d+)", basename(file)).captures
    N = parse(Int, caps[1])
    p = parse(Int, caps[2])

    # CSV einlesen
    df = CSV.read(file, DataFrame)
    y = df.y
    X = Matrix(df[:, r"^x"])

    # isample-Messungen mit @elapsed
    times = [ @elapsed matrix_regression!(X, y) for _ in 1:isample ]

    # Median und SD in µs
    median_us = median(times) * 1e6
    sd_us     = std(times)    * 1e6

    push!(results, (N=N, p=p, median_time_us=median_us, sd_time_us=sd_us))
end

# 5. Ergebnisse speichern
CSV.write(joinpath(data_path, "julia_elapsed_results.csv"), results)
