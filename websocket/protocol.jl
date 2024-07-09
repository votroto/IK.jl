using DelimitedFiles

function dlm_serialize(xs)
	buff=IOBuffer()
	writedlm(buff, xs)
	String(take!(buff))
end

function dlm_deserialize_transform(str)
	buff=IOBuffer(str)
	M = zeros(4, 4)
	M[end] = 1.0

	xs = readdlm(buff, Float64)
	M[1:3, 1:4] .= reshape(xs, 3, 4)

	@show M
	M
end