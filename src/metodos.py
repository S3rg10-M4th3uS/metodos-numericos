"""Métodos numéricos implementados com a biblioteca padrão do Python."""
import math


def finite(values):
    if not all(math.isfinite(v) for v in values):
        raise ValueError('Todos os valores devem ser finitos.')


def system(a, b):
    n = len(b)
    if not n or len(a) != n or any(len(row) != n for row in a):
        raise ValueError('A matriz deve ser quadrada e compatível com b.')
    finite([v for row in a for v in row] + list(b))
    return n


def residual(a, x, b):
    return max(abs(sum(v * w for v, w in zip(row, x)) - target)
               for row, target in zip(a, b))


def gaussian(a, b):
    """Eliminação de Gauss com pivoteamento parcial, sem alterar as entradas."""
    n = system(a, b)
    m = [list(map(float, row)) + [float(target)] for row, target in zip(a, b)]
    scale = max(abs(v) for row in a for v in row)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(m[i][k]))
        if abs(m[pivot][k]) <= 1e-14 * scale:
            raise ValueError('Sistema singular ou numericamente mal condicionado.')
        m[k], m[pivot] = m[pivot], m[k]
        for i in range(k + 1, n):
            factor = m[i][k] / m[k][k]
            m[i][k] = 0.0
            for j in range(k + 1, n + 1):
                m[i][j] -= factor * m[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][j] * x[j] for j in range(i + 1, n))) / m[i][i]
    finite(x)
    return x


def seidel(a, b, tol=1e-10, max_iter=10000, initial=None):
    n = system(a, b)
    finite([tol])
    if tol <= 0 or max_iter < 1 or any(a[i][i] == 0 for i in range(n)):
        raise ValueError('Tolerância e limite positivos; diagonal deve ser não nula.')
    x = [0.0] * n if initial is None else list(initial)
    if len(x) != n:
        raise ValueError('Vetor inicial incompatível.')
    finite(x)
    history = []
    for iteration in range(1, max_iter + 1):
        old = x.copy()
        for i in range(n):
            x[i] = (b[i] - sum(a[i][j] * x[j] for j in range(n) if i != j)) / a[i][i]
        finite(x)
        delta = max(abs(u - v) for u, v in zip(x, old))
        r = residual(a, x, b)
        history.append({'iteration': iteration, 'x': x.copy(), 'delta': delta, 'residual': r})
        if delta <= tol and r <= tol:
            return x, history
    raise ValueError(f'Gauss-Seidel não convergiu em {max_iter} iterações.')


def bridge(e, resistances, tol=1e-10):
    finite([e] + list(resistances))
    if len(resistances) != 5 or any(r <= 0 for r in resistances):
        raise ValueError('Informe cinco resistências estritamente positivas.')
    r1, r2, r3, r4, r5 = resistances
    a = [[1/r1 + 1/r4 + 1/r5, -1/r5], [-1/r5, 1/r2 + 1/r3 + 1/r5]]
    b = [e/r1, e/r2]
    (vb, vc), history = seidel(a, b, tol)
    currents = [(e-vb)/r1, (e-vc)/r2, vc/r3, vb/r4, (vb-vc)/r5]
    currents.append(currents[0] + currents[1])
    return {'matrix': a, 'rhs': b, 'voltages': [vb, vc], 'currents': currents,
            'history': history, 'residual': residual(a, [vb, vc], b)}


def least_squares(design, y):
    """Decomposição QR por Gram-Schmidt modificado com reortogonalização."""
    m, n = len(y), len(design[0])
    if m < n or len(design) != m or any(len(row) != n for row in design):
        raise ValueError('Dados insuficientes ou matriz incompatível.')
    finite([v for row in design for v in row] + list(y))
    q, r = [], [[0.0] * n for _ in range(n)]
    for j in range(n):
        v = [row[j] for row in design]
        original_norm = math.sqrt(sum(t*t for t in v))
        for _ in range(2):
            for i in range(j):
                projection = sum(t*u for t, u in zip(q[i], v))
                r[i][j] += projection
                v = [t - projection*u for t, u in zip(v, q[i])]
        r[j][j] = math.sqrt(sum(t*t for t in v))
        if r[j][j] <= 1e-13 * max(original_norm, 1):
            raise ValueError('Dados sem posto suficiente para o ajuste.')
        q.append([t/r[j][j] for t in v])
    return gaussian(r, [sum(t*u for t, u in zip(column, y)) for column in q])


def regression(x, y):
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('Informe ao menos três pares de dados.')
    finite(list(x) + list(y))
    if any(v <= 0 for v in y):
        raise ValueError('O ajuste exponencial logarítmico exige y > 0.')
    center = sum(x)/len(x)
    scale = max(abs(v-center) for v in x)
    if scale == 0:
        raise ValueError('Os valores de x devem variar.')
    t = [(v-center)/scale for v in x]
    output = {}
    for degree, name in [(1, 'reta'), (2, 'parabola')]:
        c = least_squares([[v**j for j in range(degree+1)] for v in t], y)
        pred = [sum(c[j]*v**j for j in range(degree+1)) for v in t]
        if degree == 1:
            coefficients = [c[0]-c[1]*center/scale, c[1]/scale]
        else:
            coefficients = [c[0]-c[1]*center/scale+c[2]*(center/scale)**2,
                            c[1]/scale-2*c[2]*center/scale**2, c[2]/scale**2]
        output[name] = {'coefficients': coefficients, 'predictions': pred,
                        'sse': sum((u-v)**2 for u, v in zip(y, pred))}
    alpha, beta = least_squares([[1, v] for v in t], [math.log(v) for v in y])
    a, b = math.exp(alpha-beta*center/scale), beta/scale
    pred = [math.exp(alpha+beta*v) for v in t]
    output['exponencial'] = {'coefficients': [a, b], 'predictions': pred,
                            'sse': sum((u-v)**2 for u, v in zip(y, pred)),
                            'log_sse': sum((math.log(u)-math.log(v))**2 for u, v in zip(y, pred))}
    # Ajuste adicional: minimiza o erro na escala original, sem logaritmos.
    # Para cada expoente, a amplitude ótima é calculada analiticamente.
    def profile(exponent):
        shift = max(exponent*v for v in t)
        basis = [math.exp(exponent*v-shift) for v in t]
        amplitude = sum(u*v for u,v in zip(y,basis))/sum(v*v for v in basis)
        predictions = [amplitude*v for v in basis]
        return sum((u-v)**2 for u,v in zip(y,predictions)), predictions, amplitude, shift

    middle, step = beta, 1.0
    for _ in range(60):
        left, right = middle-step, middle+step
        values = [profile(v)[0] for v in (left,middle,right)]
        if values[1] <= min(values[0], values[2]):
            break
        middle = left if values[0] < values[2] else right
        step *= 2
    else:
        raise ValueError('Não foi possível delimitar o mínimo exponencial.')
    ratio = (math.sqrt(5)-1)/2
    for _ in range(160):
        u, v = right-ratio*(right-left), left+ratio*(right-left)
        if profile(u)[0] < profile(v)[0]:
            right = v
        else:
            left = u
    exponent = (left+right)/2
    error, predictions, amplitude, shift = profile(exponent)
    b = exponent/scale
    a = math.exp(math.log(amplitude)-shift-b*center)
    output['exponencial_original'] = {'coefficients': [a,b], 'predictions': predictions, 'sse': error}
    return output


def integrate(x, y):
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('Informe ao menos três pontos para comparar os métodos.')
    finite(list(x) + list(y))
    n = len(x)-1
    h = (x[-1]-x[0])/n
    if h <= 0 or any(x[i+1] <= x[i] or not math.isclose(x[i+1]-x[i], h, rel_tol=1e-9, abs_tol=1e-12*abs(h)) for i in range(n)):
        raise ValueError('Os pontos devem ser crescentes e igualmente espaçados.')
    if n % 2:
        raise ValueError('Simpson 1/3 repetida exige número par de subintervalos.')
    return {'h': h, 'intervals': n,
            'trapezoid': h*(y[0]/2+sum(y[1:-1])+y[-1]/2),
            'simpson': h/3*(y[0]+y[-1]+4*sum(y[1:-1:2])+2*sum(y[2:-1:2]))}
