using HimalayaUI

# Package loading is the dominant startup cost (~15-20s without a sysimage).
# Exercising open_db also compiles the SQLite schema path.
mktempdir() do tmp
    HimalayaUI.open_db(joinpath(tmp, "precompile.db"))
end
