import Pkg
Pkg.add("MixedModels")


using CSV, DataFrames, MixedModels, StatsModels
using Statistics

# ---- Einstellungen ----
data_path = "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_lmm_data"
files = filter(f -> endswith(f, ".csv"), readdir(data_path; join=true))
isample = 10    # wie in deinem R-Skript

# Result-Container
results = DataFrame(model_type=String[], p=Int[], n_clusters=Int[], n_per_cl=Int[],
                    median_us=Float64[], sd_us=Float64[])

# Hilfsfunktion zur Zeitmessung
function time_fit(formula::FormulaTerm, df::DataFrame, reps::Int)
    times = Float64[]
    for _ in 1:reps
        t1 = time()
        fit(MixedModel, formula, df; REML=true)
        push!(times, time() - t1)
    end
    return times
end

# Loop
for file in files
    df = CSV.read(file, DataFrame)

    # Metadaten extrahieren
    m = match(r"lmm_(ri|rs)_p(\d+)_nc(\d+)_npc(\d+)\.csv", basename(file))
    model_label = m.captures[1]      # "ri" oder "rs"
    p_val       = parse(Int, m.captures[2])
    nc          = parse(Int, m.captures[3])
    npc         = parse(Int, m.captures[4])

    # Fixed-Teil als String
    fixed_part = join(["x$(i)" for i in 1:p_val], " + ")

    # Formel-Strings
    f1_str = "y ~ $fixed_part + (1 | cluster_id)"
    f2_str = "y ~ $fixed_part + (1 + $fixed_part | cluster_id)"

    # @formula via eval(Meta.parse(...))
    formula1 = eval(Meta.parse("@formula($f1_str)"))
    formula2 = eval(Meta.parse("@formula($f2_str)"))

    # Zeitmessung
    times_ri = time_fit(formula1, df, isample)
    times_rs = time_fit(formula2, df, isample)

    # Ergebnisse speichern
    push!(results, (
      model_type = model_label,
      p          = p_val,
      n_clusters = nc,
      n_per_cl   = npc,
      median_us  = median(times_ri)*1e6,
      sd_us      = std(times_ri)*1e6
    ))
    push!(results, (
      model_type = model_label,
      p          = p_val,
      n_clusters = nc,
      n_per_cl   = npc,
      median_us  = median(times_rs)*1e6,
      sd_us      = std(times_rs)*1e6
    ))
end

# Speichern
CSV.write("benchmark_lmm_mixedmodels.csv", results)
println("✅ Benchmark fertig und gespeichert in 'benchmark_lmm_mixedmodels.csv'")





df = CSV.read("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_lmm_data/lmm_rs_p6_nc350_npc350.csv", DataFrame)

# Modell definieren
f = @formula(y ~ x1 + x2 + x3 + x4 + x5 + x6 + (1 + x1 + x2 + x3 + x4 + x5 + x6 | cluster_id))
fit_j = fit(MixedModel, f, df)

fixef(fit_j)         # Fixe Effekte
VarCorr(fit_j)       # Random Effects
stderr(fit_j)        # Standardfehler
