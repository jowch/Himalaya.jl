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
	# peakidx, proms = peakproms(argmaxima(y), y)
	
	# smooth things out
	# k = Distributions.pdf.(Normal(0, m / 2), -(m ÷ 2):(m ÷ 2))
	# u = conv(k, y)[2:end-1]

	# ignore the initial 
	# valid_idx = remove_drop ? withoutbeamstop(u) : (1:length(u))

	u = rescale(y)

	# approximate derivatives of y
	du = savitzky_golay(m, 4, u; order = [2, 3, 4])

	# find minima in second derivative
	peakidx = let
		# find zero-crossings in third derivative	
		crossings = falses(length(u))
		crossings[1:end-1] .= diff(sign.(du[:, 2])) .!= 0

		# check if minima or maxima with fourth derivative
		# valid_idx.start .+ 
		findall(crossings .& (du[:, 3] .> 0)) .- 1
	end

	proms = let
		_, ps = peakproms(peakidx, u)
		_, qs = peakproms(peakidx .+ 1, u)
		vec(maximum(hcat(ps, qs); dims = 2))
	end

	# # estimate the second derivative of the function using Savitzky-Golay filter
	# d²u = savitzky_golay(m, n, u; order = 2)

	# # find indices of the maximum second derivative values
	# d²_idx, d²_proms = peakproms(argmaxima(-d²u), -d²u)

	if θ == 0
		# compute a θ value using the prominences of the d²u peaks
		θ = findtheta(proms)
	end

	# compute the peak prominences of peaks on the `u` scale
	θ_idx = proms .> θ
	peakidx[θ_idx], proms[θ_idx]

	# peakidx, proms
end

function findtheta(proms)
	# find a prominence threshold
	qs = [quantile(proms, p) for p = 0.75:0.01:1]

	# sanity check
	if maximum(qs) < 10 * median(qs)
		return Inf
	end

	# approximate derivatives
	dq = savitzky_golay(5, 4, qs; order = [2, 3, 4])

	inflection_idx = let
		# find zero-crossings in third derivative	
		crossings = falses(length(qs))
		crossings[1:end-1] .= diff(sign.(dq[:, 2])) .!= 0

		# check if minima or maxima with fourth derivative
		candidates = findall(crossings .& (dq[:, 3] .> 0)) .- 1

		candidates[argmax(dq[candidates, 1])]
	end

	qs[inflection_idx]
end

function withoutbeamstop(y, dy)
	findfirst(diff(sign.(dy)) .> 0):length(y)
end

function withoutbeamstop(y)
	withoutbeamstop(y, savitzky_golay(5, 1, y; order = 1))
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

