using HimalayaUI

# Package loading is the dominant startup cost (~15-20s without a sysimage).
# Exercise the SQLite schema path, the empty-args usage branch, and a real
# ArgParse dispatch (config list is side-effect-free).
mktempdir() do tmp
    HimalayaUI.open_db(joinpath(tmp, "precompile.db"))
end
HimalayaUI.main(String[])
HimalayaUI.main(["config", "list"])
