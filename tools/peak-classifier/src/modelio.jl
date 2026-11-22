function save_model(path::AbstractString, m::PeakModel)
    mkpath(dirname(abspath(path)))
    serialize(path, (version = 2, julia = string(VERSION), model = m))
    path
end

# `deserialize` runs code from the file: load only models you trained yourself.
function load_model(path::AbstractString)
    d = try
        deserialize(path)
    catch err
        error("cannot read model file $path (made by an older version of the tool?); retrain it: $(sprint(showerror, err))")
    end
    d.version == 2 || error("unsupported model file version $(d.version); retrain it")
    d.model.feature_names in (FEATURE_NAMES, filter(!in(DETECTOR_FEATURES), FEATURE_NAMES)) ||
        error("model was trained with a different feature set; retrain it")
    d.model
end
