### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ de5acd79-e83d-496d-84ab-32cce86eb135
begin
	using Pkg

	Pkg.activate("..")
	
	using Statistics
	using CSV, DataFrames, DelimitedFiles
	using Himalaya
	using Plots
	using Plots.Measures
	using LaTeXStrings
	using HypertextLiteral
	using JLD2
	using NaturalSort

	default(
		fontfamily = "Computer Modern"
	)
end;

# ╔═╡ af9d4c5c-d92a-11eb-0d27-cf0e1f21eca0
md"""
# Analyzing Sample Exposures
"""

# ╔═╡ a369fa4f-fe77-4187-8810-ed1e1f7cb1b6
DATA_DIR = "../data"

# ╔═╡ d9031353-aa48-41ed-a37b-2251c2da489f
md"""
## Under the Hood
"""

# ╔═╡ c2b64840-5750-4f2c-b914-b51567dee1ab
begin
	exposures = begin
		exposures = DataFrame(CSV.File(joinpath(DATA_DIR, "selected.csv")))
		filter(exposures) do exposure
			exposure.Use == 1
		end
	end
	
	image_dir = joinpath(DATA_DIR, "exposures")
end;

# ╔═╡ add38be0-135a-4974-b7fc-5ab1bad1a2c4
@bind target_sample @htl("""
Sample:
<select>
	$((@htl("<option value=\"$(id)\">$(id)</option>")
	   for id in sort(unique(exposures[!, :Sample]), lt = natural)))
</select>
""")

# ╔═╡ 881db234-a2b4-45dc-83b6-f494769e2756
begin
	paths = readdir(joinpath(DATA_DIR, "traces"), join=true)
	name_to_path = Dict(filename => only(filter(contains(filename), paths))
	                	for filename in exposures[!, :Filenames])
	
	function load_trace(fname)
		readdlm(name_to_path[fname], ' ', Float64, '\n')
	end
end;

# ╔═╡ ef262db0-9a2d-477d-997e-40c7333cb8c8
begin
	sample = target_sample
	fnames = exposures[exposures[!, :Sample] .== sample, :Filenames]
	traces = cat([load_trace(fname) for fname in fnames]...; dims=3)
end;

# ╔═╡ 5de1ef0c-9d87-4177-b280-6d48e75b511b
@bind target_trace @htl("""
Trace to index:
<select>
	$((@htl("<option value=\"$(id)\">$(id)</option>")
	   for id in fnames))
</select>
""")

# ╔═╡ f3bb8b27-8efb-4e1d-93b7-167b639d7333
md"## Recording Assignments"

# ╔═╡ 20fb39f0-44f5-429b-9d51-c5a8c576fa5e
begin
	assignments_path = joinpath(DATA_DIR, "assignments")
	
	if !isdir(assignments_path)
		mkdir(assignments_path)
	end

	filename = split(sample, ' ') |> first
	assignments_file = joinpath(assignments_path, "$(filename).jld2")
	
	indexing = try
		# attempt to load an existing indexing
		@info "loading $(filename) from file"
		jldopen(assignments_file, "r") do f
			f["indexing"]
		end
	catch e
		@info e
		initial = Dict(
			"sample" => sample,
			"traces" => fnames,
			# params are initialized by the default values of the input fields above
			"params" => (m = 5, n = 3, minprom = 0.001, downsample = 1),
			"assignments" => Dict(trace => [] for trace in fnames)
		)
		
		jldopen(assignments_file, "w") do f
			f["indexing"] = initial
		end
		
		initial
	end
end;

# ╔═╡ c41131e0-e6f4-4118-a3f3-541d52c45f02
@bind index_params @htl("""
<div id="param-form">
	<div id="container">
		<div class="input-box">
			<label for="m">m</label>
			<input id='m' type='number' value='$(indexing["params"][:m])' />
		</div>

		<div class="input-box">
			<label for="n">n</label>
			<input id='n' type='number' value='$(indexing["params"][:n])' />
		</div>

		<div class="input-box">
			<label for="minprom">minprom</label>
			<input id='minprom' type='number' value='$(indexing["params"][:minprom])' min='0'/>
		</div>

		<div class="input-box">
			<label for="downsample">downsample</label>
			<input id='downsample' type='number' value='$(indexing["params"][:downsample])' min='1' />
		</div>
	</div>

	<style>
		#param-form > span {
			padding: 0 1rem 0 0;
		}
		
		#container {
			display: flex;
			flex-wrap: wrap;
		}
		
		.input-box {
			width: calc(50% - 10px);
			margin: 5px;
		}

		.input-box > label {
			text-align: right;
			display: inline-block;
			width: 30%;
		}

		.input-box > input {
			width: 60%;
		}
	</style>

	<script>
		let parent = currentScript.parentElement
		let inputs = [...parent.getElementsByTagName("input")]

		let mapping = Object.fromEntries(inputs.map((x, i) => [x.id, i]))
	
		parent.value = inputs.map(input => input.value)

		inputs.forEach(input => 
			input.addEventListener('input', event => {
				parent.value[mapping[event.target.id]] = event.target.value

				parent.dispatchEvent(new CustomEvent("input"))
				event.preventDefault()
			}))
	</script>
</div>
""")

# ╔═╡ a5ad3313-baf9-40d3-b564-215be063e73d
begin
	target_trace_idx = fnames .== target_trace
	m, n, downsample = parse.(Int, index_params[[1, 2, 4]])
	minprom = parse.(Float64, index_params[3])
	
	peak_idx, proms, d2y, d2y_idx = findpeaks(
		log.(traces[:, 2, target_trace_idx]), m, n;
		minprom = minprom, downsample = downsample
	)

	peak_xs = traces[peak_idx, 1, target_trace_idx]
end;

# ╔═╡ 92b99f2c-14d6-4f46-9fa0-accbbdd07588
peak_xs

# ╔═╡ 3fe83d68-5ee4-49f1-8088-0950e3c74cc2
indices = begin
	raw_indices = indexpeaks(peak_xs; tol=0.0025)
	filtered_indices = filter(raw_indices) do index
		index.basis < 0.15
	end
	raw_indices[sortperm(score.(filtered_indices), rev = true)]
end;

# ╔═╡ b3d8d33a-499a-4fb8-8763-4e9221e44676
indices

# ╔═╡ 1b86d4d9-1fa1-471c-9d25-3fa149e252de
let
	p = plot(
		traces[:, 1, :], traces[:, 2, :], ribbon = traces[:, 3, :],
		yscale=:log10, label="I(q)"
	)
	if !isempty(peak_idx)
		vline!(p, traces[peak_idx, 1, target_trace_idx], c="gray", alpha = 0.5)
	end
	
    q = plot(d2y)

	if !isempty(d2y_idx)
    	vline!(q, d2y_idx, c="gray", alpha = 0.5)
	end
	
	r = histogram(log10.(proms), bins = 100, xticks = (-3:0.25:-1, round.(10 .^(-3:0.25:-1); digits = 4)))
    plot(p, q, r, legend = false, layout = (3, 1), size = (800, 900))
end

# ╔═╡ c09238cf-0299-4d32-b1b2-eac61a714253
log10.(proms)

# ╔═╡ 4ceb5294-3ed8-48f0-93d4-f4186874e2fc
index_info = DataFrame(Dict(
	"Index Phase" => map(x -> x.phase, indices),
	"Index Basis" => map(x -> x.basis, indices),
	"Score" => score.(indices),
	"Fit" => fit.(indices),
	"checked" => [
		index ∈ indexing["assignments"][target_trace]
		for index in indices
	]
));

# ╔═╡ ef795376-67b1-40f7-8a76-d495dc4abb7c
let
	index_check(i, name, checked) = @htl("""
		<div class="toggle-box">
			<div>
				<input
					type="checkbox" id="toggle-$i" value="$name" class="toggle-check"
					$(checked ? "checked" : nothing)
				/>
				<label class="toggle-label" for="toggle-$i">$name</label>
			</div>
		</div>
	""")

	name(row) = @htl("""
		<span>
		$(row["Index Phase"]) @ 
			<span class='name-basis'>
				$(round(row["Index Basis"]; digits = 3))
			</span>
		(
			<span class='name-score'>
				$(round(row["Score"]; digits = 3))
			</span>
		; R² = 
			<span class='name-fit'>
				$(round(row["Fit"][2]; digits = 3))
			</span>
		)
		</span>
	""")
	
	@bind target_toggles @htl("""
	<div id="trace-toggle">
		<div id="container">
			$((index_check(i, name(row), row["checked"])
			   for (i, row) in enumerate(eachrow(index_info))))
	
			<button id='clear' type='button'>Reset</button>
		</div>
	
		<style>
			#trace-toggle > span {
				padding: 0 1rem 0 0;
			}
			
			#container {
				display: flex;
				flex-wrap: wrap;
			}
			
			.toggle-box {
				width: calc(50% - 10px);
				margin: 5px;
			}
	
			.toggle-box:hover {
				background-color: #eee
			}

			.toggle-label {
				width: calc(100% - 50px);
				display: inline-block;
				margin: 0;
				padding: 0.20em;
			}

			.name-basis {
				font-weight: bold;
			}

			.name-score {
				color: #1e90ff;
			}

			.name-fit {
				color: #ff8d1e;
			}

			#container > button {
				margin: 20px 0 10px;
				width: 100%;
				padding: 0.5em 0;
			}
		</style>

		<script>
			let parent = currentScript.parentElement
			let toggles = [...parent.getElementsByTagName("input")]
				.sort((a, b) => a.id > b.id)

			let mapping = Object.fromEntries(toggles.map((x, i) => [x.id, i]))
		
			parent.value = toggles.map(toggle => toggle.checked)

			toggles.forEach(toggle => 
				toggle.addEventListener('change', event => {
					parent.value[mapping[event.target.id]] = event.target.checked

					parent.dispatchEvent(new CustomEvent("input"))
					event.preventDefault()
				}))
	
			parent.getElementsByTagName("button")[0]
				.addEventListener('click', event => {
					parent.value = toggles.map(_ => false)

					toggles.forEach(toggle => toggle.checked = false)
	
					parent.dispatchEvent(new CustomEvent("input"))
					event.preventDefault()
				})
		</script>
	</div>
	""")
end

# ╔═╡ fdf8fd82-9a29-44a7-b873-98fb112ccc12
let
	p = plot(
		traces[:, 1, :], traces[:, 2, :], ribbon=traces[:, 3, :], yscale=:log10,
		labels = reshape(fnames, 1, length(fnames)),
		size = (800, 500), xlabel = "q", ylabel = "I(q)", margin = 5mm,
		ylim = (minimum(traces[:, 2, :]) * 0.8, maximum(traces[:, 2, :] * 1.5))
	)

	if !isempty(peak_idx)
		if !any(target_toggles .> 0)
			scatter!(p,
				traces[peak_idx, 1, target_trace_idx],
				traces[peak_idx, 2, target_trace_idx] .* 1.2,
				# maximum(traces[peak_idx, 2, :], dims = 2) .* 1.2,
				marker = :dtriangle, markersize = 3, alpha = 0.75,
				markerstrokewidth = 0, color = "black", label = "Peaks"
			)
		else
			assigned_peak_idx = Set()
			
			for i in findall(target_toggles)				
				label = replace("""
				\\mathrm{$(index_info[i, "Index Phase"])},\\;
				d = $(round(index_info[i, "Fit"][1]; digits = 2)) \\;{\\AA}, 
				R^2 = $(round(index_info[i, "Fit"][2]; digits = 3))
				""", "\n" => "") |> latexstring

				peaks = collect(skipmissing(indices[i].peaks))
				
				local_peak_idx = findall(
					∈(peaks),
					skipmissing(traces[:, 1, target_trace_idx] |> vec)
				)
				
				peak_intensities = maximum(traces[local_peak_idx, 2, :], dims = 2)
				# traces[peak_idx, 2, target_trace_idx]

				scatter!(p,
					peaks,
					peak_intensities .* 1.2,
					marker = :dtriangle, markersize = 3, alpha = 0.5,
					markerstrokewidth = 0.5, label = label
				)
				
				annotate!([
					(q, qI * 1.45, (
						round(q; digits = 3), "Computer Modern", 10,
						:black, :center
					))
					for (q, qI) in zip(peaks, peak_intensities)
				])

				push!(assigned_peak_idx, local_peak_idx...)
			end
			
			unassigned_peak_idx = setdiff(peak_idx, assigned_peak_idx)
			if !isempty(unassigned_peak_idx)
				scatter!(p,
					traces[unassigned_peak_idx, 1, target_trace_idx],
					traces[unassigned_peak_idx, 2, target_trace_idx] .* 1.2,
					# maximum(traces[peak_idx, 2, :], dims = 2) .* 1.2,
					marker = :dtriangle, markersize = 3, alpha = 0.25,
					markerstrokewidth = 0, color = "black", label = "Unassigned Peak(s)"
				)
			end
		end
	end
	plot(p, legend = :topright, foreground_color_legend = nothing)
end

# ╔═╡ a173370f-aacc-43f2-8b0a-2dc815ba87b0
begin
	# update param records
	indexing["params"] = (; n, m, minprom, downsample)
	jldopen(assignments_file, "w",) do f
		f["indexing"] = indexing
	end
end

# ╔═╡ 033fc67c-8793-45bb-b965-067b56258d8a
let
	# update assignments if changed
	if any(target_toggles)
		selected_indices = indices[target_toggles]
		for index in selected_indices
			indexing["assignments"][target_trace] = selected_indices
		end
	else
		indexing["assignments"][target_trace] = []
	end
	
	jldopen(assignments_file, "w",) do f
		f["indexing"] = indexing
	end
end;

# ╔═╡ Cell order:
# ╟─af9d4c5c-d92a-11eb-0d27-cf0e1f21eca0
# ╟─a369fa4f-fe77-4187-8810-ed1e1f7cb1b6
# ╟─add38be0-135a-4974-b7fc-5ab1bad1a2c4
# ╟─5de1ef0c-9d87-4177-b280-6d48e75b511b
# ╟─c41131e0-e6f4-4118-a3f3-541d52c45f02
# ╟─fdf8fd82-9a29-44a7-b873-98fb112ccc12
# ╟─ef795376-67b1-40f7-8a76-d495dc4abb7c
# ╠═92b99f2c-14d6-4f46-9fa0-accbbdd07588
# ╠═b3d8d33a-499a-4fb8-8763-4e9221e44676
# ╟─d9031353-aa48-41ed-a37b-2251c2da489f
# ╠═de5acd79-e83d-496d-84ab-32cce86eb135
# ╠═c2b64840-5750-4f2c-b914-b51567dee1ab
# ╠═881db234-a2b4-45dc-83b6-f494769e2756
# ╠═ef262db0-9a2d-477d-997e-40c7333cb8c8
# ╠═a5ad3313-baf9-40d3-b564-215be063e73d
# ╠═3fe83d68-5ee4-49f1-8088-0950e3c74cc2
# ╠═1b86d4d9-1fa1-471c-9d25-3fa149e252de
# ╠═c09238cf-0299-4d32-b1b2-eac61a714253
# ╠═4ceb5294-3ed8-48f0-93d4-f4186874e2fc
# ╟─f3bb8b27-8efb-4e1d-93b7-167b639d7333
# ╠═20fb39f0-44f5-429b-9d51-c5a8c576fa5e
# ╠═a173370f-aacc-43f2-8b0a-2dc815ba87b0
# ╠═033fc67c-8793-45bb-b965-067b56258d8a
