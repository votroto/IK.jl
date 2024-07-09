using HTTP.WebSockets

desired = """
0.6418551170035527	-0.3919805344677002	0.6590700033947623	0.06719912719304007
-0.3729819081066398	0.5913721269480282	0.7149569942969264	0.6562296440221194
-0.6700048544611344	-0.70471999266899	0.23337357804844314	0.8664256519532297
0.0	0.0	0.0	1.0
"""

WebSockets.open("ws://127.0.0.1:8081") do sock
	send(sock, desired)
	println("waiting...")
	angles = receive(sock)
	println(angles)
end
