import numpy as np, json
import matplotlib.pyplot as plt

#

PRECISION = 8
norma = lambda x: abs(x)

EPS = 10 ** -PRECISION

# pentru un x dat, aprox f(x) folosind:
# * Newton progresiv si Aitken
# * cmm patrate, Horner

def read_data(file_name = None):
	global decl_f
	if file_name is None:
		f = decl_f
		x0, xn, n = float(input('x0=')), float(input('xn=')), int(input('n='))
		if x0 > xn:
			x0, xn = xn, x0
		h = (xn - x0)/(n-1)
		xs = [x0+i*h for i in range(n-1)] + [xn]
		fxs = [f(x) for x in xs]
		return xs, fxs, f

	data = json.load(open(file_name, "rt"))
	if data['method'] == 'given':
		return data['xs'], data['fxs'], lambda x: 0
	elif data['method'] == 'compute':
		exec(data['f'])
		x0, xn, n = data['x0'], data['xn'], data['n']
		if x0 > xn:
			x0, xn = xn, x0
		h = (xn - x0)/(n-1)
		xs = [x0+i*h for i in range(n-1)] + [xn]
		fxs = [f(x) for x in xs]
		return xs, fxs, f
	else:
		raise Exception("Method not good")


def Newton_prog(xs, fxs):
	n = len(xs)
	h = (xs[-1] - xs[0])/(len(xs)-1)

	# Schema Aitken
	ait = [[] for _ in range(n)]
	ait[0] = [fxs[i+1]-fxs[i] for i in range(n-1)]
	for i in range(1, n-1):
		ait[i] = [ait[i-1][j+1]-ait[i-1][j] for j in range(n-1-i)]
	# for i in range(n-1):
	# 	print(ait[i])

	def Ln(x):
		t = (x - xs[0])/h
		ss = [t]
		for k in range(2, n):
			ss.append(ss[-1]*(t-k+1)/k)

		# print(prev, fxs[0])
		l = fxs[0]
		l += sum([ait[k][0]*ss[k] for k in range(n-1)])
		return l
	return Ln

def cmm_patrate(xs, fxs):
	global M
	return np.linalg.solve()


def Horner(polinom, x):
	d = 0
	for i in range(len(polinom)):
		d = polinom[i] + d * x
	return d

def plott(xs, fxs, f_aprox, actual_f):
	plt.plot(xs, [f_aprox(x) for x in xs], linewidth=5, color='orange', label='f_aprox')
	plt.plot(xs, [actual_f(x) for x in xs], linewidth=1, color='blue', label='actual')
	# plt.plot(xs, fxs, linewidth=1, color='blue', label='actual')
	plt.legend()
	plt.show()

def decl_f(x):
	return x ** 4 - 12 * (x ** 3) + 30 * (x ** 2) + 12

M = 5
xs, fxs, actual_f = read_data('t6.json')
print("xs, fxs:", xs, fxs)
ll = Newton_prog(xs, fxs)
x_ = xs[0] + (xs[1] - xs[0])/2
# print(ll(x_), actual_f(x_))

print(f"{ll(x_)=}\n{norma(ll(x_) - actual_f(x_))=}")
plott(xs, fxs, ll, actual_f)