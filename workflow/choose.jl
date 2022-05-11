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

# ╔═╡ 2915e3c4-000d-4d3f-9fe6-2414e6e9377e
begin
	using Pkg

	Pkg.activate("..")
	
	using CSV, DataFrames, DelimitedFiles, Plots
	using PlutoUI, HypertextLiteral
	
	TableOfContents()
end

# ╔═╡ 9ed9ce8b-1e68-4984-80c4-ddba2f76dc42
pwd()

# ╔═╡ b5b65743-8b43-4272-8ed2-cd3157fadd37
@bind data_dir html"<input type='text' />"

# ╔═╡ 39ef1558-4fa3-454d-af99-9d626b7ed366
@bind samples_file html"<input type='text' />"

# ╔═╡ 3cc6c21c-e6a9-4ad1-b17a-61eefe8ebb7a
begin
	exposures = DataFrame(CSV.File(joinpath(data_dir, samples_file)))
	
	if !("Use" in names(exposures))
		exposures[!, :Use] .= 1
	end
end;

# ╔═╡ 5d774c81-5906-43dd-afae-46ff386ed548
@bind target_sample @htl("""
<select>
	$((@htl("<option value=\"$(id)\">$(id)</option>")
	   for id in unique(exposures[!, :Sample])))
</select>
""")

# ╔═╡ fdca0ea6-d211-471b-9a68-4333c1be395f
begin
	using Base64, Images, TiffImages
	using FixedPointNumbers, Colors
	
	image_dir = joinpath(data_dir, "exposures")
	
	image_paths = filter(contains(".tif"), readdir(image_dir, join=true))
	img_path_map = Dict(filename => only(filter(contains(filename), image_paths))
						for filename in exposures[!, :Filenames])
	
	function load_image(fname; exposure_dir=image_dir, width=nothing)
		img = reinterpret(Gray{N0f32}, TiffImages.load(img_path_map[fname]))
		
		a = log.(Array{Float32}(img) .+ eps())
		a_eq = adjust_histogram((a .- minimum(a)) ./ (maximum(a) - minimum(a)),
						    	Images.Equalization())
		
		if !isnothing(width)
			new_size = trunc.(Int, size(a_eq) .* (width / size(a_eq, 2)))
			a_eq = imresize(a_eq, new_size)			
		end
		
		Gray.(convert.(N0f32, clamp.(a_eq, 0, 1)))
	end
	
	function image_to_blob(image)
		stringmime("image/png", image)
	end
end;

# ╔═╡ 6689cd27-83b5-4f0d-a980-ebc8042c98f3
md"""
# Exposure Choice

This tool is designed to help you review the exposures that you took for each sample and decide which ones you want to keep for further analysis.
"""

# ╔═╡ 77984666-b203-4a73-985e-b410b13cd9ef
md"""
To start, please enter the path of the data directory in the box below. This folder should contain a samples file, a exposures file, and a folder called `tot_files`.  The path that you enter can be an absolute path or a relative path. For reference here is the current directory that this notebook is running in:
"""

# ╔═╡ 246e66bc-6fcc-43c5-afd3-e748aae1f6f5
md"""
## Selecting the Index Files

Enter the path (relative to the `data` folder you entered above) to a CSV file containing all of the exposures you took at the beamline. This file should mimic the structure of our lab notebook. The column names that are expected are `Sample`, `Time (s)`, and `Filename`.
"""

# ╔═╡ 3e0aeb28-6216-4889-bff2-1f6bd51a5892
md"""
## Reviewing Exposures

You can choose which sample you want to review using the dropdown selector below. Once selected, the graph below will update to displace the tot traces for all (selected) exposures of this sample. You can use the check boxes below to toggle each trace on and off. Whether or not the traces are selected will be marked in a column called `Use` in the exposures index CSV file.
"""

# ╔═╡ c4288e70-2d4e-4178-9409-4bc9b3c547a5
md"## Under the Hood"

# ╔═╡ ac4445ce-268d-40a9-bf29-1822d7fb886c
begin
	paths = readdir(joinpath(data_dir, "traces"), join = true)
	path_map = Dict(filename => only(filter(contains(filename), paths))
					for filename in exposures[!, :Filenames])
end;

# ╔═╡ 6b243905-4e92-47d3-a28a-2b3766ecf03b
target_exposures = exposures[exposures[!, :Sample] .== target_sample, :Filenames];

# ╔═╡ 2b844eea-2ede-4e1f-ac94-66c9f11c5324
begin
	targets = exposures[in(target_exposures).(exposures.Filenames),
						[:Filenames, :Use]];
	
	targets[!, :Preview] = [image_to_blob(load_image(fname, width=100))
							for fname in targets[!, :Filenames]]
end;

# ╔═╡ a5a1f6c7-7fd7-45b5-b765-34c18d79511c
begin
	exposure_check(name, checked, blob) = @htl("""
		<div class="toggle-box">
			<img src="data:image/png;base64,$(blob)" />
			<div>
				<input
					type="checkbox" id="$name" value="$name"
					$(checked ? "checked" : nothing)
				/>
				<label for="$name">$name</label>
			</div>
		</div>
	""")
	
	@bind target_toggles @htl("""
	<div id="trace-toggle">
		<div id="container">
			$((exposure_check(name, checked == 1, blob)
			   for (name, checked, blob) in eachrow(targets)))	
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
				width: calc(25% - 10px);
				margin: 5px;
			}
		
			.toggle-box > img {
				width: 100%;
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
		</script>
	</div>
	""")
end

# ╔═╡ 5b5ad09c-ecbc-4e4c-8b08-eb41e70600ab
begin
	target_tots = cat([readdlm(path_map[path], ' ', Float64, '\n')
			   	   	   for path in target_exposures]...; dims=3)
	
	plot(
		target_tots[:, 1, target_toggles], target_tots[:, 2, target_toggles],
		ribbon = target_tots[:, 3, target_toggles],
		yaxis = :log, title = "Selected Traces for #$target_sample",
		labels = reshape(
			target_exposures[target_toggles], 1,
			length(target_exposures[target_toggles])
		),
		dpi = 300
	)
end

# ╔═╡ f251b6c9-8587-4412-a4cc-c23f57f17a71
begin
	for (name, value) in zip(target_exposures, target_toggles)
		exposures[exposures[!, :Filenames] .== name, :Use] .= value
	end
	
	CSV.write(joinpath(data_dir, "selected.csv"), exposures)
end;

# ╔═╡ Cell order:
# ╟─6689cd27-83b5-4f0d-a980-ebc8042c98f3
# ╟─77984666-b203-4a73-985e-b410b13cd9ef
# ╟─9ed9ce8b-1e68-4984-80c4-ddba2f76dc42
# ╟─b5b65743-8b43-4272-8ed2-cd3157fadd37
# ╟─246e66bc-6fcc-43c5-afd3-e748aae1f6f5
# ╟─39ef1558-4fa3-454d-af99-9d626b7ed366
# ╟─3e0aeb28-6216-4889-bff2-1f6bd51a5892
# ╟─5d774c81-5906-43dd-afae-46ff386ed548
# ╟─5b5ad09c-ecbc-4e4c-8b08-eb41e70600ab
# ╟─a5a1f6c7-7fd7-45b5-b765-34c18d79511c
# ╟─c4288e70-2d4e-4178-9409-4bc9b3c547a5
# ╠═2915e3c4-000d-4d3f-9fe6-2414e6e9377e
# ╠═3cc6c21c-e6a9-4ad1-b17a-61eefe8ebb7a
# ╠═ac4445ce-268d-40a9-bf29-1822d7fb886c
# ╠═6b243905-4e92-47d3-a28a-2b3766ecf03b
# ╠═fdca0ea6-d211-471b-9a68-4333c1be395f
# ╠═2b844eea-2ede-4e1f-ac94-66c9f11c5324
# ╠═f251b6c9-8587-4412-a4cc-c23f57f17a71
