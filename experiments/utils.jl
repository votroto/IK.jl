function check_obj(sol, w, init) 
    sum(2 .* w .* (1 .- cos.(sol) .* cos.(init) .- sin.(sol) .* sin.(init)))
end