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
function findpeaks(y; θ = 0, m = 5, n = 3)
	u = y .* rescale(y)
	
	# estimate the second derivative of the function using Savitzky-Golay filter
	d²u = let
		d²u = savitzky_golay(m, n, u; nd = 2)
		d²u[d²u .>= 0] .= 0
		d²u
	end

	# find indices of the maximum second derivative values
	d²_idx, d²_proms = peakproms(argmaxima(-d²u), -d²u)

	if θ == 0
		# compute a θ value using the prominences of the d²u peaks
		θ = findtheta(d²_proms)
	end

	# compute the peak prominences of peaks on the `u` scale
	peakproms(d²_idx[d²_proms .>= θ], u)
end

function findtheta(proms)
	# find a prominence threshold
	qs = [quantile(proms, p) for p = 0.75:0.01:1]

	# sanity check
	if maximum(qs) < 10 * median(qs)
		return Inf
	end

	d²q = let
		d²q = savitzky_golay(5, 3, log10.(qs); nd = 2)
		d²q[d²q .< 0] .= 0
		d²q
	end

	qs[argmax(d²q)]
end

"""
Computes the Savitzky-Golay filtering of `y`.
"""
function savitzky_golay(m, n, y; nd = 0)
	num_y = length(y)
	z = -m:m
	J = zeros(2m + 1, n + 1)

	for i = 0:n
		@inbounds J[:, i+1] .= z .^ i
	end

	# The convolution term matrix 
	C = J' \ I(n .+ 1)[:, nd .+ 1] # = inv(J' * J) * J' = pinv(J)
	Y = zeros(num_y, length(nd))

	for i in 1:num_y
		if i <= m
			window_indices = abs.(z .+ i) .+ 1
		elseif i > num_y - m
			window_indices = -abs.(z .+ i .- num_y) .+ num_y
		else
			window_indices = z .+ i
		end

		for j in eachindex(nd)
			@inbounds Y[i, j] = C[:, j]' * y[window_indices]
		end
	end

	length(nd) == 1 ? Y[:, 1] : Y
end



