import math, random

m = -1
u = 0.1

while 1 + u != 1:
	m -= 1
	u *= 0.1
u*=10
m += 1


print(f"Ex1:\n{u=}, {m=}")

x, y, z = 1.0, u/10, u/10
print(f"\nEx2:")
print(f"{x, y, z=}")
print(f"{((x + y) + z != x + (y + z))=}")

m = -1
u = 0.1

x, y, z = 1.4587634893721988298, 9.43983289434897328449821, 2349.39454795868969865895679856656565
print(f"{x, y, z=}")
print(f"{((x * y) * z != x * (y * z))=}")

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
			return (10395 * a - 1260 * a ** 3 + 21 * a ** 5) / (10395 - 4725 * a ** 2 + 210 * a ** 4 - a ** 6)
		case 7:
			return (135135 * a - 17325 * a ** 3 + 378 * a ** 5 - a ** 7) / (135135 - 62370 * a ** 2 + 3150 * a ** 4 - 28 * a ** 6)
		case 8:
			return (2027025 * a - 270270 * a ** 3 + 6930 * a ** 5 - 36 * a ** 7) / (2027025 - 945945 * a ** 2 + 51975 * a ** 4 - 630 * a ** 6 + a ** 8)
		case 9:
			return (34459425 * a - 4729725 * a ** 3 + 135135 * a ** 5 - 990 * a ** 7 + a ** 9) / (34459425 - 16216200 * a ** 2 + 945945 * a ** 4 - 13860 * a ** 6 + 45 * a ** 8)
def S(n, a):
	return T(n, a) / math.sqrt(1 + (T(n, a) ** 2))
def C(n, a):
	return 1 / math.sqrt(1 + (T(n, a) ** 2))

def test_10_000(f, fverif):
	nrs = [random.uniform(-math.pi/2, math.pi/2) for _ in range(10_000)]
	rez = [() for _ in range(10_000)]
	for i in range(4, 10):
		for j in range(len(nrs)):
			err = abs(f(i, j) - fverif(j))
			if rez[j] == () or rez[j][0] > err:
				rez[j] = (err, i)

	r = [list(map(lambda x: x[0], filter(lambda x: x[1] == i, rez))) for i in range(4, 10)]

	l = [(sum(r[i - 4])/len(r[i - 4]), i) for i in range(4, 10)]
	l.sort(key=lambda x: x[0])

	for i in l:
		print(f"{f.__name__}({i[1]}, a) are eroarea medie {i[0]}")
	print("=====")

	return l[0][1] # returneaza cea mai buna valoare pentru n

print("\nEx3:")

best_tan = test_10_000(T, math.tan)
best_sin = test_10_000(S, math.sin)
best_cos = test_10_000(C, math.cos)
print(f"{best_tan=}\n{best_sin=}\n{best_cos=}\n=====")

nr = random.uniform(-math.pi/2, math.pi/2)

print(f"{nr=}")
print(f"T({best_tan}, nr) = {T(best_tan, nr):<12} | {math.tan(nr)=} (err={abs(math.tan(nr) - T(best_tan, nr))})")
print(f"S({best_sin}, nr) = {S(best_sin, nr):<12} | {math.sin(nr)=} (err={abs(math.sin(nr) - S(best_sin, nr))})")
print(f"C({best_cos}, nr) = {C(best_cos, nr):<12} | {math.cos(nr)=} (err={abs(math.cos(nr) - C(best_cos, nr))})")