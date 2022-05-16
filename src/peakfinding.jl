using Peaks
using LinearAlgebra

"""
	findpeaks(trace, m, n)

Identifies peaks in the trace using the second derivative of the Savitzky-Golay.
"""
function findpeaks(trace; m = 5, n = 3, top = 15)
	# estimate the second derivative of the function using Savitzky-Golay filter
	d2y = let
		d2y = savitzky_golay(m, n, trace; nd = 2)
		d2y[d2y .>= 0] .= 0
		d2y
	end

	# peaks in trace
	idx, proms = peakproms(argmaxima(-d2y), -d2y)
	top_k = sortperm(proms; rev = true)[1:top]	

	idx[top_k], proms[top_k]
end

"""
Computes the second derivative of the Savitzky-Golay filtering of `y`.
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



