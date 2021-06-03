import numpy as np
import matplotlib.pyplot as plt
import os

# defining the parameters
r = 0.05
u = 1. + 0.15
d = 1 - 0.05
p = (1 + r - d) / (u - d)
q = (u - (1 + r)) / (u - d)
S_0 = 100
K = 100
N = 5

S = np.empty(shape=(N + 1, N + 1))
S.fill(np.nan)
S[0, 0] = S_0

for n in range(1, N + 1):
    for k in range(0, n + 1):
        S[n, k] = np.round(S_0 * np.power(d, n - k) * u ** k, 3)

time = np.argwhere((~np.isnan(S)))[:, 1].T
S_1d = S[np.logical_not(np.isnan(S))]

print('S\n', S)


APP = np.empty(shape=(N + 1, N + 1))
APP.fill(np.nan)
APP[0, 0] = 0

for i in range(1, N + 1):
    for j in range(0, N + 1):
        APP[i, j] = np.max([S_0 - S[i, j], 0])

print('APP\n', APP)
APP_1d = APP[np.logical_not(np.isnan(APP))]

APEP = np.empty(shape=(N + 1, N + 1))
APEP.fill(np.nan)
APEP[N, :] = APP[N, :]

# M[i, j], i riga, j colonne
for i in range(N, 0, -1):
    for j in range(0, i):
        APEP[i - 1, j] = (q * APEP[i, j] + p * APEP[i, j + 1]) / (1 + r)

APEP_1d = APEP[np.logical_not(np.isnan(APEP))]
print('APEP \n', APEP)

Omega = np.empty(shape=(N, 2 ** N))
Omega.fill(np.nan)

for n in range(N):
    for m in range(2 ** N):
        if (m % 2 ** (n + 1)) < 2 ** n:
            Omega[N - n - 1, m] = 0
        else:
            Omega[N - n - 1, m] = 1

print('Omega \n', Omega)
APMV = np.maximum(APP, APEP)

APMV_1d = APMV[np.logical_not(np.isnan(APMV))]
print('APMV \n', APMV)

OS = np.empty(shape=(N + 1, 2 ** N))
OS.fill(np.nan)
OS[N, :] = np.repeat(N, 2 ** N)
list(range(4, 0, -1))
prova = np.zeros(((N-1)*2**N, 3))
i = 0
for n in range(N - 1, 0, -1):
    for m in range(2 ** N):
        if APMV[n, int(Omega[0:(n), m].sum())] == np.max([APP[n, int(Omega[0:(n), m].sum())], 0]):
            OS[n, m] = n
        else:
            OS[n, m] = OS[n + 1, m]

OS[0, :] = 0 if APMV[0,0] == APP[1,1] else OS[1, :]