# Globally Optimal Inverse Kinematics

Use via a primitive websocket server
```sh
julia websocket/server.jl -i 0.0.0.0 -p 8081 data/kuka.toml
```
then you can send the first 3x4 elements of a transformation matrix.
```sh
cat data/feasible_kuka_pose.txt | websocat ws://0.0.0.0:8081
```

The server will remember the last pose and minimize the squared distance from it.