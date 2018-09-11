import numpy as np


def coef(x, y):
    '''
    This function computes Newton's coefficients to evaluate the polynomial
    that passes all the data points (xi, yi) in vector x and vector y.
    Inputs:
        x : array of data points
        y : array of f(x)
    Output:
        polyCoeff: array of coefficients used to evaluate f() based on x.
                   Note: it is NOT the polynomial's coefficient vector.
    '''
    x = np.array(x)
    y = np.array(y)
    x.astype(float)
    y.astype(float)
    n = len(x)
    polyCoeff = y[:n]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            polyCoeff[i] = float(polyCoeff[i] -
                                 polyCoeff[i - 1]) / float(x[i] - x[i - j])

    return np.array(polyCoeff)  # return an array of coefficients


def eval(polyCoeff, x, r):
    '''
    This function computes f(r) based on polyCoeff and x.
    Inputs:
        polyCoeff: see above
        x: see above. It must be the vector used to compute polyCoeff.
        r: could be a scalar or an array of scalars, to be evaluated.
    Outputs:
        value: = f(r)

    '''
    x = np.array(x)
    x.astype(float)
    n = len(polyCoeff) - 1
    value = polyCoeff[n]
    r = np.array(r)
    for i in range(n - 1, -1, -1):
        value = value * (r - x[i]) + polyCoeff[i]
    return value  # return the y_value interpolation


def example():
    x = np.array([1, 2, 3])
    y = np.array([1, 4, 9])
    r = np.array([4, 5, 6])
    polyCoeff = coef(x, y)
    value = eval(polyCoeff, x, r)
    print('data points: ', r)
    print('evaluation: ', value)
    assert (value == [16, 25, 36]).all()


example()
