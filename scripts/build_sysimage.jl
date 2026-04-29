using PackageCompiler

repo_root = dirname(@__DIR__)
out_dir   = joinpath(repo_root, "build")
out_path  = joinpath(out_dir, "himalaya.so")
workload  = joinpath(repo_root, "scripts", "precompile_workload.jl")

mkpath(out_dir)

@info "Building sysimage → $out_path (this takes several minutes)"

create_sysimage(
    ["HimalayaUI"];
    sysimage_path             = out_path,
    precompile_execution_file = workload,
    project                   = joinpath(repo_root, "packages", "HimalayaUI"),
)

write(joinpath(out_dir, "julia_version"), string(VERSION))

@info """
Sysimage built: $out_path

Usage:
  julia --project=packages/HimalayaUI --sysimage $(relpath(out_path, repo_root)) \\
        -e 'using HimalayaUI; main(ARGS)' -- serve /path/to/experiment
"""
