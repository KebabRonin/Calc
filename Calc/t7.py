import random, math
kmax = 30_000
def get_r(p: list[float]):
	aa = max(map(abs, p))
	return (abs(p[0]) + aa) / abs(p[0])

def bonus(f, f_der):
	iters=200
	EPS = 1e-5
	x = random.random()
	y = f(x)
	for i in range(iters):
		if abs(f(x)) <= EPS:
			return x
		if abs(f_der(x)) <= EPS:
			break
		if abs(f_der(x)*(f(x)-f(y))) <= EPS:
			break
		y = x - f(x)/f_der(x)
		x = x - (f(x)**2+f(y)**2)/(f_der(x)*(f(x)-f(y)))
	return x


def Horner(polinom, x):
	d = 0
	for i in range(len(polinom)):
		d = polinom[i] + d * x
	return d

sign = lambda x: -1 if x <= 0 else 1

def Mueller(polinom, x0, x1, x2):
	k = 3
	dx = 10
	while EPS <= abs(dx) <= 1e8 and k <= kmax:
		h0 = x1 - x0
		h1 = x2 - x1
		d0 = (Horner(polinom, x1) - Horner(polinom, x0)) / h0
		d1 = (Horner(polinom, x2) - Horner(polinom, x1)) / h1
		a = (d1 - d0) / (h0 + h1)
		b = a*h1 + d1
		c = Horner(polinom, x2)

		if b**2-4*a*c<0:
			break
		if abs(b+sign(b)*math.sqrt(b**2-4*a*c)) <= EPS:
			break

		dx = (2*c)/(b+sign(b)*math.sqrt(b**2-4*a*c))
		xk = x2 - dx
		k += 1
		x0, x1, x2 = x1, x2, xk
	# print("Iters:", k)
	if abs(dx) < EPS:
		# solutie
		return x0
	else:
		# divergenta
		return None

EPS = 1e-15
polinom = [1, -6, 11, -6]
polinom = [42, -55, -42, 49, -6]
polinom = [8, -38, 49, -22, 3]
polinom = [1, -6, 13, -12, 4]

f=lambda x:math.exp(x)-math.sin(x)
f_der=lambda x: math.exp(x)-math.cos(x)
nr_incercari = 30_000

print(bonus(f, f_der))
exit(0)

r = get_r(polinom)
print("R=", r)
# print(Horner(polinom, 2/3))
# exit(0)
sols = []
for _ in range(nr_incercari):
	x0, x1, x2 = random.uniform(-r, r), random.uniform(-r, r), random.uniform(-r, r)
	# print("Inceput:", x0, x1, x2)
	s = Mueller(polinom, x0, x1, x2)
	if s is not None:
		# print("Solutia:", s)
		for ss in sols:
			if abs(s - ss) <= 1e-5:
				break
		else:
			sols.append(s)
	else:
		# print("Divergenta")
		...


sols.sort()
print("Nr sols:", len(sols))
print(sols)
sols = list(map(lambda x: x.__str__() + '\n', sols))
with open("t7.out", "wt") as f:
	f.writelines(sols)