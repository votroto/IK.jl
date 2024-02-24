using LinearAlgebra

function solve_forward_kinematics(x, d, r, α)
    prod(dh_t.(x, d, α, r))
end

function random_feasible_pose(d, r, α, θl, θh)
    x = θl .+ rand(length(r)) .* (θh .- θl)
    prod(dh_t.(x, d, α, r))
end

function random_bounded_pose(xmin, xmax, ymin, ymax, zmin, zmax)
    function _rotmat(w, theta)
        bw = [0 -w[3] w[2]; w[3] 0 -w[1]; -w[2] w[1] 0]
        Matrix(I, 3, 3) + sin(theta) * bw + (1 - cos(theta)) * bw * bw
    end

    function _inner(d, r, α, θl, θh)
        R = _rotmat(normalize(rand(3)), (rand() - 0.5) * 2pi)
        x = xmin + rand() * (xmax - xmin)
        y = ymin + rand() * (ymax - ymin)
        z = zmin + rand() * (zmax - zmin)

        [
            R[1, 1] R[1, 2] R[1, 3] x
            R[2, 1] R[2, 2] R[2, 3] y
            R[3, 1] R[3, 2] R[3, 3] z
            0 0 0 1
        ]
    end
    _inner
end

function pose_error(A, B)
    pos_e = norm(A[1:3, 4] - B[1:3, 4])
    rot_e = acos(clamp(0.5 * (tr(A[1:3, 1:3] \ B[1:3, 1:3]) - 1), -1.0, 1.0))

    pos_e, rot_e
end
