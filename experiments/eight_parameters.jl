using LinearAlgebra

r = rand(8)
d = rand(8)
α = rand(8) .- 0.5
θl = -rand(8)
θh = rand(8)

w = normalize(ones(8), 1)
θ = zeros(8)


#=
r = [0.3909466035900754, 0.05334042898961133, 0.7718202347243546, 0.628982030818105, 0.8687922030938143, 0.18285393976530495, 0.5358776362220101, 0.5668543950485855]
d = [0.9325325959247484, 0.4458894136088887, 0.040835928167835434, 0.5582525626918645, 0.42230900145155825, 0.8026322846539268, 0.7264171915898002, 0.5500498447819858]
α = [0.3695885790273007, 0.48346033503441077, 0.17890793185414555, 0.0491638641181571, 0.2423971093501599, 0.291212328657674, 0.30151417676089354, -0.05881750680771436]
θl = [-0.07905983486738777, -0.20511794089659974, -0.23731514110785878, -0.5283727597868716, -0.9558468076884618, -0.49947371118058126, -0.31941685950031584, -0.03637867240821502]
θh = [0.8683352697717882, 0.862057635136576, 0.6323443159677039, 0.8268104923736423, 0.5495375630997179, 0.7661245471180052, 0.06834142500972129, 0.7433185004182264]
w = [0.35355339059327373, 0.35355339059327373, 0.35355339059327373, 0.35355339059327373, 0.35355339059327373, 0.35355339059327373, 0.35355339059327373, 0.35355339059327373]
θ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
=#


#=fast?=#
 r=[
 0.4097466815667613
 0.3028629154089517
 0.4923686842167555
 0.8299960052771359
 0.9830510767694165
 0.028025488949107502
 0.4616594461342548
 0.28164160409799666
 ]

d=[
 0.9214296616331262
 0.43446742676691386
 0.4661439157908347
 0.3773575439906235
 0.6769107477177506
 0.9325010411902442
 0.02890523372665066
 0.9283389119959244
]
α=[

  0.09712414283244386
 -0.01623887516185374
  0.3457782787898437
  0.2081874646962497
  0.4892317152763235
  0.004065525778685353
  0.3787493186751357
 -0.48663389548643365
]
θl=[
 -0.6713572204768526
 -0.1489518772928593
 -0.1764350566273406
 -0.8466384137056324
 -0.08614218362111437
 -0.12627052591977117
 -0.9687195259470519
 -0.5648650630776294
]
θh = [
 0.7348353733615461
 0.8678692309086151
 0.8674938095496425
 0.40357490814170993
 0.1995446551584913
 0.022277816071957335
 0.25463035330831163
 0.7693469365366229
]