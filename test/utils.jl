rand_between(lo, hi) = rand() * (hi - lo) + lo
rand_between(dims...; lo=-1, hi=1) = rand(dims...) .* (hi - lo) .+ lo
