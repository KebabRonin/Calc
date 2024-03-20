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
PRECISION = 6
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

def calc_b(A_init, s):
	b = numpy.zeros_like(s)
	for i in range(s.shape[0]):
		b[i] = sum([s[j] * A_init[i, j] for j in range(s.shape[0])])
	return b

def QR_decomp(A_init, b_init=None):
	if b_init is not None:
		b = copy.deepcopy(b_init)
	else:
		b = None
	Qt = numpy.identity(A_init.shape[0])
	R = copy.deepcopy(A_init)
	n = A_init.shape[0]
	for r in range(n-1):
		sigma = sum([R[i, r] ** 2 for i in range(r, n)])
		if is_zero(sigma):
			print("Matricea A e singulara")
			break
		k = math.sqrt(sigma) * (-1 if R[r, r] > 0 else 1)
		beta = sigma - k * R[r, r]

		u = numpy.zeros((n,))
		u[r] = R[r, r] - k
		for i in range(r+1, n):
			u[i] = R[i, r]

		# R = P * R
		for j in range(r+1, n):
			gamma = sum([u[i] * R[i, j] for i in range(r, n)]) / beta
			for i in range(r, n):
				R[i, j] -= gamma * u[i]
		R[r, r] = k
		for i in range(r+1, n):
			R[i, r] = 0

		# b = P * b
		if b is not None:
			gamma = sum([u[i] * b[i] for i in range(r, n)]) / beta
			for i in range(r, n):
				b[i] -= gamma * u[i]

		# Q = P * Q
		for j in range(n):
			gamma = sum([u[i] * Qt[i, j] for i in range(r, n)]) / beta
			for i in range(r, n):
				Qt[i, j] -= gamma * u[i]

	Q = numpy.transpose(Qt)
	return Q, R

def subinv(A,b):
	n = A.shape[0]
	x = numpy.zeros_like(b)
	for i in reversed(range(n)):
		x[i] = (b[i] - sum([A[i, j] * x[j] for j in range(i+1, n)])) / A[i, i]
	return x

def solve(A_init, b_init, Q, R):
	Qt = numpy.transpose(Q)
	n = A_init.shape[0]
	bQ = (Qt @ numpy.atleast_2d(b_init).T)
	for i in range(n):
		if is_zero(R[i, i]):
			print("R nu are det 0, ecuatia nu are solutie")
			return None
	x = subinv(R, bQ)
	return x

def QRinv(A_init):
	Q, R = QR_decomp(A_init)
	n = A_init.shape[0]
	for i in range(n):
		if is_zero(R[i, i]):
			print("R nu are det 0, inversa nu se poate calcula")
			return None
	A_inv = numpy.zeros_like(A_init)
	for j in range(n):
		b = numpy.array([Q[j, i] for i in range(n)])
		x_prim = subinv(R, b)
		for i in range(n):
			A_inv[i, j] = x_prim[i]
	return A_inv


def bonus(A_init):
	n = A_init.shape[0]
	for i in range(n):
		for j in range(i+1, n):
			if A_init[i, j] != A_init[j, i]:
				print("[BONUS] Matricea nu e simetrica", i, j)
				return
	mat_c = copy.deepcopy(A_init)
	mat_next = numpy.zeros_like(A_init)
	last_err = 10
	while not is_zero(last_err):#norma(mat_next - mat_c)
		Q, R = QR_decomp(mat_c)
		mat_next = R @ Q
		if norma(mat_next - mat_c) == last_err:
			print("Stuck in a loop")
			return None
		else:
			last_err = norma(mat_next - mat_c)
		print(last_err, norma(mat_next - mat_c))
		print(mat_c, "\n\n", mat_next, "\n\n")
		mat_c = mat_next
	return mat_next

def afisare(A_init, s):
	b_init = calc_b(A_init, s)
	n = s.shape[0]
	libs = numpy.linalg.qr(A_init)
	xQR = solve(A_init, b_init, libs.Q, libs.R)
	tems = QR_decomp(A_init, b_init)
	xHouseholder = solve(A_init, b_init, tems[0], tems[1])
	print("Lib:")
	print("Q=\n", libs.Q)
	print("R=\n", libs.R)
	print("xQR=\n", xQR)
	print("Tema:")
	print("Q=\n", tems[0])
	print("R=\n", tems[1])
	print('xHouseholder=\n', xHouseholder)
	err = norma(xQR - xHouseholder)
	print(f"norma(xQR - xHouseholder): {err}")

	print("----> Ex4: Erori")
	err1 = norma(A_init @ xHouseholder - numpy.atleast_2d(b_init).T)
	err2 = norma(A_init @ xQR - numpy.atleast_2d(b_init).T)
	err11 = norma(xHouseholder - numpy.atleast_2d(s).T) / norma(s)
	err21 = norma(xQR - numpy.atleast_2d(s).T) / norma(s)
	print("norma(A * xH  - b):         ", err1)
	print("norma(A * xQR - b):         ", err2)
	print("norma(A * xH  - s)/norma(s):", err11)
	print("norma(A * xQR - s)/norma(s):", err21)
	print("inv:\n", QRinv(A_init))
	print("norma(invH - invNumpy):     ", norma(QRinv(A_init) - numpy.linalg.inv(A_init)))

A_init, s = read_data(IN_FILE)
# A_init, s = rand_data(5)
afisare(A_init, s)
# bonus_ex = numpy.matrix([[1, 6, 2], [6, 3, 4], [2, 4, 10]])
# print(bonus(bonus_ex))