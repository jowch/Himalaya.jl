# Usage (from tools/peak-classifier):
#   julia --threads=auto --project=. scripts/peak_classifier.jl cv    DB [--group-by experiment|sample]
#   julia --threads=auto --project=. scripts/peak_classifier.jl train DB --out models/lab.jls
#   julia --threads=auto --project=. scripts/peak_classifier.jl predict MODEL TRACE.dat
# Shared options for cv and train: --threshold 0.5 --nms 8 --trees 200 --no-form-factor --no-detector-features
using PeakClassifier, HimalayaUI

function parse_opts(args)
    pos = String[]; opt = Dict{String, String}()
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--no-form-factor", "--no-detector-features")
            opt[a[3:end]] = "1"
        elseif startswith(a, "--")
            i += 1; i <= length(args) || error("$a needs a value")
            opt[a[3:end]] = args[i]
        else
            push!(pos, a)
        end
        i += 1
    end
    pos, opt
end

function main(args)
    isempty(args) && error("usage: peak_classifier.jl cv|train|predict ...")
    cmd, rest = args[1], args[2:end]
    pos, opt = parse_opts(rest)
    kw = (threshold = parse(Float64, get(opt, "threshold", "0.5")),
          nms = parse(Int, get(opt, "nms", "8")),
          n_trees = parse(Int, get(opt, "trees", "200")),
          detector_features = !haskey(opt, "no-detector-features"))
    if cmd == "predict"
        length(pos) == 2 || error("predict needs MODEL and TRACE.dat")
        m = load_model(pos[1])
        q, I, σ = HimalayaUI.load_dat(pos[2])
        foreach(println, predict_peaks(m, q, I, σ))
        return
    end
    length(pos) == 1 || error("$cmd needs a database path")
    t0 = time()
    exps = load_labeled_exposures(pos[1]; group_by = Symbol(get(opt, "group-by", "experiment")),
                                  include_form_factor = !haskey(opt, "no-form-factor"))
    println("loaded $(length(exps)) exposures in $(round(time() - t0; digits = 1)) s")
    ds = build_dataset(exps)
    println("candidates: $(length(ds.label)) rows; ",
            join(["$(n)=$(count(==(l), ds.label))" for (n, l) in
                  ("positive" => 1, "negative" => 0, "weak_negative" => -1, "dropped" => -2)], ", "))
    if cmd == "cv"
        PeakClassifier.print_cv(cross_validate(exps, ds; kw...))
    elseif cmd == "train"
        haskey(opt, "out") || error("train needs --out PATH")
        m = train(ds; kw...)
        save_model(opt["out"], m)
        println("saved $(opt["out"]) ($(m.n_train.positive) positive, $(m.n_train.negative) negative)")
    else
        error("unknown command $cmd")
    end
    println("total $(round(time() - t0; digits = 1)) s")
end

main(ARGS)
