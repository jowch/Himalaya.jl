module DataFramesExt

using HimalayaDB, DataFrames, Tables

HimalayaDB.dataframe(rows) = DataFrame(Tables.columntable(rows))

end # module
