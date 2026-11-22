module HimalayaDB

using SQLite, DBInterface, Tables, SparseArrays
import Himalaya

export connect,
    experiments, setups, samples, exposures,
    curated_peaks, index_candidates, confirmed_indices,
    phase_call, sample_labels, subtraction_background, bridged_samples,
    reconstruct_index, load_trace, dataframe

"""
    dataframe(rows) -> DataFrame

Convert Tables.jl rows to a DataFrame. Method provided by the DataFrames
weakdep extension — call `using DataFrames` to activate it.
"""
function dataframe end

include("connect.jl")
include("queries.jl")
include("sample_annotations.jl")
include("reconstruct.jl")
include("trace.jl")

end # module
