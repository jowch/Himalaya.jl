using Lazy

function int2miller(n)
    bits = bitstring(n)
    first_nonzero = findfirst(!x -> x == '0', bits)
    padding = 3 - ((length(bits) - first_nonzero + 1) % 3)

    minimal_leading = bits[(first_nonzero - padding):end]

    triples = [parse.(Int8, split(minimal_leading[i:i+2], ""))
               for i in 1:3:length(minimal_leading)]

    reduce((a, b) -> a .+ b, triples)
end

miller_indices = @lazy distinct(map(x -> int2miller(x), Lazy.range()));
