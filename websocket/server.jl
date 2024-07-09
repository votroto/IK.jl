include("../src/IK.jl")

include("parse_manipulator.jl")
include("protocol.jl")

using LinearAlgebra
using HTTP.WebSockets
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--port", "-p"
            help = "Websocket listening port"
			arg_type = Int
            default = 8081
        "--ip", "-i"
			help = "Websocket listening ip"
            default = "127.0.0.1"
        "dh_file"
            help = "TOML file of the DH parameters"
            required = true
    end

    parse_args(s)
end

function main(dh_file, ip, port)
	params = parse_manipulator_file(dh_file)
	d, r, α, θl, θh = params
	n = length(r)

	w = ones(n) ./ n
	θ = θl .+ rand(length(d)) .* (θh .- θl)

	WebSockets.listen(ip, port) do sock
		for msg in sock
			M = dlm_deserialize_transform(msg)
			θ, = IK.solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w)
			send(sock, dlm_serialize(θ))
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	arg = parse_commandline()
    main(arg["dh_file"], arg["ip"], arg["port"])
end
