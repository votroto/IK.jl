using TOML

function parse_manipulator_data(data)
    parsed =  TOML.parse(data)

    r = parsed["length"]
    d = parsed["offset"]
    a = parsed["angle"]

    @assert length(r) == length(d) == length(a)
    nr = length(r)

    lb = get(parsed, "lower", fill(-π, nr))
    ub = get(parsed, "upper", fill(+π, nr))

    d, r, a, lb, ub
end

function parse_manipulator_file(file)

        data = read(file, String)
        parse_manipulator_data(data)

end