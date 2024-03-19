import numpy, math, sys, os, copy

# 1. Calc b
# 2. Descompunerea QR
# 3. Rezolvarea sistemului
# 4. Erori
# 5. Inversa A
# 6. Random init
# Bonus: limita A^(k)
# Restrictii: nr mic matrici, inmultirea cu partea superioara R, fara if-uri


# Init: precizia si fisierul date ca argumente la linia de comanda
PRECISION = 10
IN_FILE = "t3.in"
norma = lambda x: numpy.linalg.norm(x, ord=2)

for arg in sys.argv[1:]:
	try:
		PRECISION = int(arg)
	except ValueError:
		if os.path.isfile(arg):
			IN_FILE = arg

EPS = 10 ** -PRECISION

def read_data(file_name: str):
	with open(file_name, "rt") as f:
		data = f.read()

	data = list(map(lambda a: float(a), data.split()))

	n = int(data[0])
	data = data[1:]

	if len(data) != n ** 2 + n:
		raise Exception("Invalid data format. Format: first number is n(matrix dimensions), next n * n are matrix values, next n are s(vector) values")

	s = numpy.asarray(data[n ** 2:])
	A_init = numpy.matrix([[x for x in data[i*n:i*n+n]] for i in range(n)])
	return A_init, s

def rand_data(n):
	A_init = numpy.random.rand(n, n)
	s = numpy.random.rand(n)
	return A_init, s


def is_zero(x: float) -> bool:
	return math.fabs(x) < EPS

def sols(A_init, s):
	n = s.shape[0]
	print(f"\n----> Ex3: Solutii:")
	# x = solve(A_init, s)
	# x_bon = solve_bonus(A_init, s)
	x_lib = numpy.linalg.solve(A_init, s)
	A_inv_lib = numpy.linalg.inv(A_init)
	# err1 = norma(A_init * x.reshape((n, 1)) - s)
	err11 = norma(A_init * x_lib.reshape((n, 1)) - s)
	# err2 = norma(x - x_lib)
	# err3 = norma(x - A_inv_lib * s.reshape((n, 1)))
	# print(f"Solutie calc : {x} ca numpy?: {is_zero(norma(x - x_lib))}")
	# print(f"Solutie bonus: {x_bon} ca numpy?: {is_zero(norma(x_bon - x_lib))}")
	print(f"Solutie numpy: {x_lib}")
	print(f"")
	# print(f"{'norma(A_init * x - s)':<35}: {err1:<25}; {(err1  < 1e-8)=}")
	print(f"{'norma(A_init * x_lib - s)':<35}: {err11:<25}; {(err11 < 1e-8)=}")
	# print(f"{'norma(x - x_lib)':<35}: {err2:<25};")
	# print(f"{'norma(x - A_inv_lib * s)':<35}: {err3:<25};")
	print(f"Inversa A^-1 numpy:\n{A_inv_lib}")

def afisare(A_init, s):
	sols(A_init, s)

def calc_b(A_init, s):
	b = numpy.zeros_like(s)
	for i in range(s.shape[0]):
		b[i] = sum([s[j] * A_init[i, j] for j in range(s.shape[0])])
	return b

def QR_decomp(A_init, b_init):
	b = copy.deepcopy(b_init)
	Q = numpy.identity(A_init.shape[0])
	R = copy.deepcopy(A_init)
	n = A_init.shape[0]
	for r in range(n - 1):
		sigma = sum([A_init[j, r] ** 2 for j in range(r, n)])
		if is_zero(sigma):
			print("Matricea A e singulara")
			break
		k = math.sqrt(sigma) * (1 if A_init[r, r] >= 0 else -1)
		beta = sigma + k * A_init[r, r]
		u = numpy.zeros((n,))
		u[r] = A_init[r, r] - k
		for i in range(r+1, n):
			u[i] = A_init[i, r]
		f = 1 / beta
		# build Pr
		P = numpy.zeros_like(A_init)
		for i in range(r):
			P[i, i] = 1
		for i in range(r, n):
			for j in range(r, n):
				A_init[i, j] = -f * u[i] * u[j]
			A_init[i, i] += 1
		R = P * R
		b = P * b.reshape(3, 1)
		Q = P * Q
	return Q, R

A_init, s = read_data(IN_FILE)
print(QR_decomp(A_init, calc_b(A_init, s)))
# afisare(A_init, s)