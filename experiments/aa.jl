# number of permanent
p = 5
# number of non-permanent
n = 6

# number of swings of permanent
#   every OTHER perm. agree * at least 2 non-perm. agree
sp = binomial(p - 1, p - 1) * sum(binomial(n, i) for i in 2:n)
# number of swings of non-permanent
#   all perm. agree * exactly one of the OTHER non-perm. agrees
sn = binomial(p, p) * binomial(n - 1, 1)

# normalized-banzhaf of permanent
bp = sp / (sp * p + sn * n)
# normalized-banzhaf of non-permanent
bn = sn / (sp * p + sn * n)