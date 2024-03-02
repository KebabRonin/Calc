import math, random

m = -1
u = 0.1

while 1 + u != 1:
	m -= 1
	u *= 0.1
u*=10
m += 1


print(f"Ex1: {u=}, {m=}")

x, y, z = 1.0, u/10, u/10
print(f"Ex2: {(x + y) + z != x + (y + z)=}")

def T(t, a):
	match t:
		case 1:
			return a
		case 2:
			return (3 * a) / (3 - a ** 2)
		case 3:
			return (15 * a - a ** 3) / (15 - 6 * a ** 2)
		case 4:
			return (105 * a - 10 * a ** 3) / (105 - 45 * a ** 2 + a ** 4)
		case 5:
			return (945 * a - 105 * a ** 3 + a ** 5) / (945 - 420 * a ** 2 + 15 * a ** 4)
		case 6:
			return (10395 - 1260 * a ** 3 + 21 * a ** 5) / (10395 - 4725 * a ** 2 + 210 * a ** 4 - a ** 6)
		case 7:
			return (135135 * a - 17325 * a ** 3 + 378 * a ** 5 - a ** 7) / (135135 - 62370 * a ** 2 + 3150 * a ** 4 - 28 * a ** 6)
		case 8:
			return (2027025 * a - 270270 * a ** 3 + 6930 * a ** 5 - 36 * a ** 7) / (2027025 - 945945 * a ** 2 + 51975 * a ** 4 - 630 * a ** 6 + a ** 8)
		case 9:
			return (34459425 * a - 4729725 * a ** 3 + 135135 * a ** 5 - 990 * a ** 7 + a ** 9) / (34459425 - 16216200 * a ** 2 + 945945 * a ** 4 - 13860 * a ** 6 + 45 * a ** 8)


nrs = [random.uniform(-math.pi/2, math.pi/2) for _ in range(10_000)]
rez = [() for _ in range(10_000)]


for i in range(4, 10):
	for j in range(len(nrs)):
		err = abs(T(i, j) - math.tan(j))
		if rez[j] == () or rez[j][0] > err:
			rez[j] = (err, i)

r = [list(map(lambda x: x[0], filter(lambda x: x[1] == i, rez))) for i in range(4, 10)]

l = [(sum(r[i])/len(r[i]), i) for i in range(4, 10)]
l.sort(key=lambda x: x[1])

for i in l:
	print(f"T({i}, a) are eroarea medie {l[i][0]}")