using Peaks
using LinearAlgebra

"""
	findpeaks(trace, m, n)

Identifies peaks in the trace using the second derivative of the Savitzky-Golay
"""
function findpeaks(trace, m, n; minprom = 1, downsample = 1)
	# estimate the second derivative of the function using Savitzky-Golay filter
	downsample_indices = 1:downsample:length(trace)
	nd2y = let 
		d2y = savitzky_golay(m, n, trace[downsample_indices]; nd = 2)
		nd2y = zeros(size(d2y))
		nd2y[findall(d2y .< 0)] .= d2y[findall(d2y .< 0)]

		nd2y
	end

	# peaks in trace
	best_indices, proms = peakproms(argminima(nd2y), nd2y; minprom = minprom)
	
	downsample_indices[best_indices], proms, nd2y, best_indices
end

"""
Computes the second derivative of the Savitzky-Golay filtering of `y`.
"""
function savitzky_golay(m, n, y; nd = 0)
	num_y = length(y)
	z = -m:m
	J = zeros(2m + 1, n + 1)

	for i = 0:n
		J[:, i+1] .= z .^ i
	end

	# The convolution term matrix 
	C = J' \ I(n .+ 1)[:, nd .+ 1] # = inv(J' * J) * J' = pinv(J)

	Y = zeros(num_y, length(nd))

	@inbounds for i in 1:num_y
		if i <= m
			window_indices = abs.(z .+ i) .+ 1
		elseif i > num_y - m
			window_indices = -abs.(z .+ i .- num_y) .+ num_y
		else
			window_indices = z .+ i
		end

		@inbounds for j in eachindex(nd)
			Y[i, j] = C[:, j]' * y[window_indices]
		end
	end

	length(nd) == 1 ? Y[:, 1] : Y
end



