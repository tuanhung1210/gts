import numpy as np
from numpy import linalg as LA
from scipy.linalg import svd

# Input A

data = np.genfromtxt('input.txt', delimiter=' ')
A = np.array(data)

m = len(A)
n = len(A[0])
AT = A.T

# Right Singular Vectors
if n > m:
    res = np.dot(AT, A)
    w, v = LA.eig(res)
    i = 0
    while (i < len(w) - 1):
        j = i + 1
        while j < len(w):
            if w[i] < w[j]:
                w[[i, j]] = w[[j, i]]
                v[:, [i, j]] = v[:, [j, i]]
            j = j + 1
        i = i + 1

# Left Singular Vectors
else:
    res = np.dot(A, AT)
    w, u = LA.eig(res)
    i = 0
    while i < len(w) - 1:
        j = i + 1
        while j < len(w):
            if w[i] < w[j]:
                w[[i, j]] = w[[j, i]]
                u[:, [i, j]] = u[:, [j, i]]
            j = j + 1
        i = i + 1

# Middle Singular Vectors
s = []
for i in range(n):
    temp = [0] * m
    s.append(temp)
for i in range(min(n, m)):
    s[i][i] = np.sqrt(w[i])
s = np.array(s)

# Cal the other Singular Vectors
if n > m:
    u = []
    for i in range(m):
        t = []
        for j in range(n):
            t.append(v[j][i])
        t = (1 / s[i][i]) * (A.dot(t))
        u.append(t)
    u = np.array(u)
    u = u.T
else:
    v = []
    for i in range(n):
        t = []
        for j in range(m):
            t.append(u[j][i])
        t = (1 / s[i][i]) * (AT.dot(t))
        v.append(t)
    v = np.array(v)
    v = v.T

# Print Ans
print("A=")
print(A)
print("U = ")
print(u)
print("S = ")
print(s.T)
print("V = ")
print(v)



