using Distributions
# using DSP: conv

"""
	findpeaks(y; θ, m, n)

Identifies peaks in curve `y` using the Savitzky-Golay second derivative.

The values of `y` are first weighted by the relative magnitude of their values.
This is to emphasize the higher value peaks and to minimize the lower ones as
these are more likely to be from noise. The second derivative of this weighted
curve 

Most of the peaks identified in this way will come from noise in the data and
will have near-zero prominence. Real peaks will have prominences higher than
those of the noisy peaks.

We can compute a threshold prominence value ``θ``, where peaks with prominence
greater than ``θ`` are considered real. This threshold can be found by looking
at the precentile-prominence curve, which will have "jumps" at higher
percentiles as real peaks have higher prominence than noise. The method for
computing this threshold is a heuristic and was developed by studying
percentile-prominence curves from many traces.

See also `Peaks.peakproms`.
"""
function findpeaks(y; θ = 0, m = 5)
	# rescale the y values
	u = y ./ maximum(y)

	# find maxima in the rescaled curve
	peaks = Himalaya.Peaks.findmaxima(u)
	peaks = Himalaya.Peaks.peakproms(peaks)

	# approximate 2nd derivative of y
	d2u = savitzky_golay(m, 4, u; order = 2)
	d2u ./= maximum(abs.(d2u))

	# find minima in the 2nd derivative
	d2peaks = let
		(; indices, heights) = Himalaya.Peaks.findminima(d2u)
		valid_idx = heights .<= -0.005

		(; indices = indices[valid_idx], proms = heights[valid_idx])
	end
	d2peaks = Himalaya.Peaks.peakproms(d2peaks)

	# find indices of d2peaks that are the same or one away from a peak in peaks
	index_distances = abs.(peaks.indices .- d2peaks.indices')
	overlapping_indices = Tuple.(findall(index_distances .<= 1))

	# use the peaks of the original trace for the overlapping peaks
	overlapping_peaks = peaks.indices[first.(overlapping_indices)]
	# use the larger prominence for the overlapping peaks
	overlapping_proms = map(overlapping_indices) do (i, j)
		max(peaks.proms[i], -d2peaks.proms[j])
	end

	peaks, proms = let upeaks = peaks
		peak_idx = setdiff(1:length(upeaks.indices), first.(overlapping_indices))
		d2peak_idx = setdiff(1:length(d2peaks.indices), last.(overlapping_indices))

		peaks = vcat(upeaks.indices[peak_idx], d2peaks.indices[d2peak_idx], overlapping_peaks)
		proms = vcat(upeaks.proms[peak_idx], -d2peaks.proms[d2peak_idx], overlapping_proms)
		
		sort_idx = sortperm(peaks)

		peaks[sort_idx], proms[sort_idx]
	end

	# compute the peak prominences of peaks on the `u` scale
	θ_idx = proms .≥ θ
	peaks[θ_idx], proms[θ_idx]
end

"""
	savitzky_golay(m, n, y; order = 0)
	savitzky_golay(m, n, y; order = [1, 2, 3])

Computes the Savitzky-Golay filtering of `y`.
"""
function savitzky_golay(m, n, y; order = 0)
	# @assert n >= order "`n` must be at least `order`"

	num_y = length(y)
	z = -m:m
	J = zeros(2m + 1, n + 1)

	for i = 0:n
		@inbounds J[:, i+1] .= z .^ i
	end

	# The convolution term matrix 
	C = J' \ I(n .+ 1)[:, order .+ 1] # = inv(J' * J) * J' = pinv(J)
	Y = zeros(num_y, length(order))

	for i in 1:num_y
		if i <= m
			window_indices = abs.(z .+ i) .+ 1
		elseif i > num_y - m
			window_indices = -abs.(z .+ i .- num_y) .+ num_y
		else
			window_indices = z .+ i
		end

		for j in eachindex(order)
			@inbounds Y[i, j] = C[:, j]' * y[window_indices]
		end
	end

	length(order) == 1 ? Y[:, 1] : Y
end

