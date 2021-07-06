using Peaks

function find_peaks(trace; minprom=5)
	# peaks in trace
	idx, peaks = peakprom(trace, 2, minprom=minprom)
	
	idx, peaks
	
# 	# peaks in derivative
# 	dtrace = mapwindow(mean, diff(trace), 9, "reflect")
# 	didx, dpeaks = peakprom(dtrace .* 5, 5, minprom=minprom)
	
# 	@show didx
	
# 	both_idx = vcat(idx, didx)
# 	order = sortperm(both_idx)
	
# 	both_idx[order], vcat(peaks, trace[didx])[order]
end
