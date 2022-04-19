v = -2

def chi2(x):
    return (x[0] - 1)**2 + (x[1] + v)**2 + phi(x[2])

def phi(t):
        return (t+1)**2