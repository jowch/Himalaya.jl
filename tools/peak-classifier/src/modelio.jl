const MODEL_FILE_VERSION = 3

function save_model(path::AbstractString, m::TrainedModel)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON.json(io, Dict{String, Any}("file_version" => MODEL_FILE_VERSION, "model" => to_dict(m)))
    end
    path
end

function load_model(path::AbstractString)
    d = try
        JSON.parse(read(path, String))
    catch err
        error("cannot read model file $path (made by an older version of the tool?); retrain it: $(sprint(showerror, err))")
    end
    (d isa AbstractDict && haskey(d, :file_version) && d[:file_version] == MODEL_FILE_VERSION) ||
        error("unsupported model file version; retrain it")
    m = from_dict(TrainedModel, d[:model])
    m.feature_version == feature_version(m.recipe.features) ||
        error("model was trained with features \"$(m.feature_version)\" but this build has \"$(feature_version(m.recipe.features))\"; retrain it")
    m
end
