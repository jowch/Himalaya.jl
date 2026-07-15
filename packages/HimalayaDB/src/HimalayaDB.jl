module HimalayaDB

using SQLite, DBInterface, Tables, SparseArrays
import Himalaya

export connect,
    experiments, samples, exposures,
    curated_peaks, index_candidates, confirmed_indices,
    reconstruct_index, load_trace, dataframe

"""
    dataframe(rows) -> DataFrame

Convert Tables.jl rows to a DataFrame. Method provided by the DataFrames
weakdep extension — call `using DataFrames` to activate it.
"""
function dataframe end

include("connect.jl")
include("queries.jl")
include("reconstruct.jl")
include("trace.jl")

end # module
