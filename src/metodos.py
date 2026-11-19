"""Métodos numéricos com exibição detalhada dos cálculos."""
import math


# --- Funções auxiliares de formatação para exibir matrizes e vetores ---

def _format_vector(vec, decimals=6):
    """Formata uma lista de números para exibição limpa no terminal."""
    items = ', '.join(f'{val:.{decimals}f}' for val in vec)
    return f'[{items}]'


def _format_augmented_matrix(matrix, decimals=4):
    """Formata a matriz aumentada [A | b] alinhada por colunas."""
    lines = []
    for row in matrix:
        coefs = '  '.join(f'{val:10.{decimals}f}' for val in row[:-1])
        rhs = f'{row[-1]:10.{decimals}f}'
        lines.append(f'  | {coefs} | {rhs} |')
    return '\n'.join(lines)


# --- Validações e cálculo de resíduo ---

def finite(values):
    """Garante que nenhum valor na lista seja NaN ou Infinito."""
    if not all(math.isfinite(val) for val in values):
        raise ValueError('Todos os valores devem ser finitos.')


def system(a, b):
    """Valida se a matriz A é quadrada, compatível com b e possui apenas números finitos."""
    n = len(b)
    if n == 0 or len(a) != n:
        raise ValueError('A matriz deve ser quadrada e compatível com b.')

    for row in a:
        if len(row) != n:
            raise ValueError('A matriz deve ser quadrada e compatível com b.')
        finite(row)

    finite(b)
    return n


def residual(a, x, b):
    """Calcula o maior erro (resíduo máximo) entre A*x e b."""
    max_error = 0.0
    for row, target in zip(a, b):
        approx = sum(coef * val for coef, val in zip(row, x))
        error = abs(approx - target)
        if error > max_error:
            max_error = error
    return max_error


# --- 1. Eliminação de Gauss ---

# Define a função que recebe a matriz de coeficientes (a), o vetor de demandas (b)
# e o parâmetro verbose (que ativa ou desativa a impressão detalhada na tela)
def gaussian(a, b, verbose=True):
    """Resolve Ax = b por Eliminação de Gauss com pivoteamento parcial."""
    # Confere se 'a' é quadrada, compatível com 'b' e sem valores inválidos, devolvendo o tamanho n
    n = system(a, b)

    # Cria uma cópia da matriz aumentada [A | b] convertendo tudo para decimal (float) sem alterar os dados originais
    m = [list(map(float, row)) + [float(target)] for row, target in zip(a, b)]
    # Descobre o maior número (em módulo) de toda a matriz para usar como régua de escala
    scale = max(abs(val) for row in a for val in row)

    # Se verbose for True, imprime a matriz aumentada inicial formatada em colunas antes de começar as contas
    if verbose:
        print('\n--- ELIMINAÇÃO DE GAUSS ---')
        print('Matriz aumentada inicial [A | b]:')
        print(_format_augmented_matrix(m))

    # Etapa 1: Escalonamento — percorre cada coluna k da diagonal principal (de 0 até n - 1)
    for k in range(n):
        # Começa supondo que o maior pivô da coluna k já está na própria linha diagonal k
        pivot = k
        # Percorre todas as linhas i abaixo da diagonal k procurando o maior valor em módulo
        for i in range(k + 1, n):
            if abs(m[i][k]) > abs(m[pivot][k]):
                # Atualiza o índice pivot caso encontre um número maior nessa coluna
                pivot = i

        # Verifica se o maior pivô disponível é zero ou praticamente zero em relação à escala da matriz
        if abs(m[pivot][k]) <= 1e-14 * scale:
            # Interrompe o cálculo e avisa que o sistema não possui solução única
            raise ValueError('Sistema singular ou numericamente mal condicionado.')

        # Exibe qual foi o pivô escolhido e em qual linha ele está
        if verbose:
            print(f'\n[Etapa {k + 1}] Coluna {k + 1} — Pivô escolhido: {m[pivot][k]:.6f} (linha {pivot + 1})')

        # Se o maior pivô estiver em outra linha, troca a linha k com a linha pivot
        # Isso evita divisão por zero e reduz erros de arredondamento no computador
        if pivot != k:
            m[k], m[pivot] = m[pivot], m[k]
            if verbose:
                print(f'  -> Trocando linha {k + 1} com linha {pivot + 1}.')

        # Percorre cada linha i abaixo do pivô para zerar o elemento da coluna k
        for i in range(k + 1, n):
            # Calcula o multiplicador de linha dividindo o elemento a eliminar pelo pivô
            factor = m[i][k] / m[k][k]
            # Zera diretamente a posição abaixo do pivô para não sobrar resíduo decimal mínimo
            m[i][k] = 0.0
            # Atualiza as colunas restantes à direita (incluindo o vetor b): Li = Li - factor * Lk
            for j in range(k + 1, n + 1):
                m[i][j] -= factor * m[k][j]

            # Mostra a operação matemática realizada nesta linha
            if verbose:
                print(f'  -> L{i + 1} = L{i + 1} - ({factor:.6f}) * L{k + 1}')

        # Exibe como a matriz ficou após zerar toda a coluna k
        if verbose and k < n - 1:
            print('  Matriz após eliminação:')
            print(_format_augmented_matrix(m))

    # Exibe a matriz triangular superior pronta e inicia a etapa final
    if verbose:
        print('\nMatriz triangular superior final:')
        print(_format_augmented_matrix(m))
        print('\n--- SUBSTITUIÇÃO RETROATIVA ---')

    # Etapa 2: Substituição retroativa — cria a lista x de tamanho n para guardar as respostas
    x = [0.0] * n
    # Percorre as equações de baixo para cima (da última linha n - 1 até a primeira 0)
    for i in range(n - 1, -1, -1):
        # Soma os termos à direita da diagonal que já possuem valores de x conhecidos
        known_sum = sum(m[i][j] * x[j] for j in range(i + 1, n))
        # Isola a incógnita x[i]: (lado direito - soma dos conhecidos) / coeficiente da diagonal
        x[i] = (m[i][n] - known_sum) / m[i][i]

        # Mostra no terminal a conta feita para descobrir esta incógnita
        if verbose:
            print(f'  x[{i + 1}] = ({m[i][n]:.6f} - {known_sum:.6f}) / {m[i][i]:.6f} = {x[i]:.8f}')

    # Faz uma verificação final de segurança para garantir que nenhum resultado virou Infinito ou NaN
    finite(x)
    # Retorna a lista x com os volumes calculados de cada mina
    return x

# --- 2. Gauss-Seidel e Ponte de Wheatstone ---

def seidel(a, b, tol=1e-10, max_iter=10000, initial=None, verbose=True):
    """Resolve Ax = b por tentativas sucessivas (método iterativo de Gauss-Seidel)."""
    n = system(a, b)
    finite([tol])

    if tol <= 0 or max_iter < 1 or any(a[i][i] == 0 for i in range(n)):
        raise ValueError('Tolerância e limite positivos; diagonal deve ser não nula.')

    x = [0.0] * n if initial is None else list(initial)
    if len(x) != n:
        raise ValueError('Vetor inicial incompatível.')
    finite(x)

    if verbose:
        print('\n--- MÉTODO DE GAUSS-SEIDEL ---')
        print(f'Chute inicial x(0) = {_format_vector(x)} | Tolerância = {tol:.2e}')

    history = []

    for iteration in range(1, max_iter + 1):
        old_x = x.copy()

        for i in range(n):
            other_terms = sum(a[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (b[i] - other_terms) / a[i][i]

        finite(x)

        delta = max(abs(new_val - old_val) for new_val, old_val in zip(x, old_x))
        err = residual(a, x, b)
        history.append({
            'iteration': iteration,
            'x': x.copy(),
            'delta': delta,
            'residual': err,
        })

        if verbose:
            converged = (delta <= tol and err <= tol)
            if iteration <= 10 or iteration % 50 == 0 or converged:
                print(
                    f'  Iteração {iteration:3d}: x = {_format_vector(x, 8)} '
                    f'| Passo (delta) = {delta:.3e} | Resíduo = {err:.3e}'
                )
            elif iteration == 11:
                print('  ... (calculando iterações intermediárias) ...')

        if delta <= tol and err <= tol:
            if verbose:
                print(f'Convergência atingida em {iteration} iterações!')
            return x, history

    raise ValueError(f'Gauss-Seidel não convergiu em {max_iter} iterações.')


def bridge(e, resistances, tol=1e-10, verbose=True):
    """Calcula tensões nodais e correntes de uma Ponte de Wheatstone."""
    finite([e] + list(resistances))
    if len(resistances) != 5 or any(r <= 0 for r in resistances):
        raise ValueError('Informe cinco resistências estritamente positivas.')

    r1, r2, r3, r4, r5 = resistances
    g1, g2, g3, g4, g5 = [1 / r for r in resistances]

    a = [
        [g1 + g4 + g5, -g5],
        [-g5, g2 + g3 + g5],
    ]
    b = [e * g1, e * g2]

    if verbose:
        print('\n--- PONTE DE WHEATSTONE (ANÁLISE NODAL) ---')
        print(f'1. Condutâncias (G = 1/R): G1={g1:.6f}, G2={g2:.6f}, G3={g3:.6f}, G4={g4:.6f}, G5={g5:.6f} S')
        print('2. Equações nodais (Leis de Kirchhoff):')
        print(f'   Nó B: ({a[0][0]:.6f})*VB + ({a[0][1]:.6f})*VC = {b[0]:.6f}')
        print(f'   Nó C: ({a[1][0]:.6f})*VB + ({a[1][1]:.6f})*VC = {b[1]:.6f}')

    (vb, vc), history = seidel(a, b, tol, verbose=verbose)

    i1 = (e - vb) / r1
    i2 = (e - vc) / r2
    i3 = vc / r3
    i4 = vb / r4
    i5 = (vb - vc) / r5
    i6 = i1 + i2
    currents = [i1, i2, i3, i4, i5, i6]

    if verbose:
        print('\n3. Cálculo das correntes (Lei de Ohm: I = ΔV / R):')
        print(f'   I1 = ({e} - {vb:.6f}) / {r1} = {i1:.8f} A')
        print(f'   I2 = ({e} - {vc:.6f}) / {r2} = {i2:.8f} A')
        print(f'   I3 = {vc:.6f} / {r3} = {i3:.8f} A')
        print(f'   I4 = {vb:.6f} / {r4} = {i4:.8f} A')
        print(f'   I5 = ({vb:.6f} - {vc:.6f}) / {r5} = {i5:.8f} A')
        print(f'   I6 (Total) = I1 + I2 = {i1:.8f} + {i2:.8f} = {i6:.8f} A')

    return {
        'matrix': a,
        'rhs': b,
        'voltages': [vb, vc],
        'currents': currents,
        'history': history,
        'residual': residual(a, [vb, vc], b),
    }


# --- 3. Mínimos Quadrados e Regressão ---

def least_squares(design, y, verbose=False):
    """Resolve mínimos quadrados via decomposição QR (Gram-Schmidt modificado)."""
    m = len(y)
    if not design or m < len(design[0]):
        raise ValueError('Dados insuficientes ou matriz incompatível.')

    n = len(design[0])
    if len(design) != m or any(len(row) != n for row in design):
        raise ValueError('Dados insuficientes ou matriz incompatível.')

    for row in design:
        finite(row)
    finite(y)

    q = []
    r = [[0.0] * n for _ in range(n)]

    for j in range(n):
        col = [row[j] for row in design]
        original_norm = math.sqrt(sum(val * val for val in col))

        for _ in range(2):
            for i in range(j):
                proj = sum(q_val * c_val for q_val, c_val in zip(q[i], col))
                r[i][j] += proj
                col = [c_val - proj * q_val for c_val, q_val in zip(col, q[i])]

        norm = math.sqrt(sum(val * val for val in col))
        r[j][j] = norm

        if norm <= 1e-13 * max(original_norm, 1):
            raise ValueError('Dados sem posto suficiente para o ajuste.')

        q.append([val / norm for val in col])

    qt_y = [sum(q_val * y_val for q_val, y_val in zip(q_col, y)) for q_col in q]

    if verbose:
        print('  [QR] Matriz triangular R:')
        print(_format_augmented_matrix([row + [rhs] for row, rhs in zip(r, qt_y)]))

    return gaussian(r, qt_y, verbose=False)


def regression(x, y, verbose=True):
    """Ajusta reta, parábola, exponencial logarítmica e exponencial direta aos pontos."""
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('Informe ao menos três pares de dados.')

    finite(list(x) + list(y))
    if any(val <= 0 for val in y):
        raise ValueError('O ajuste exponencial logarítmico exige y > 0.')

    center = sum(x) / len(x)
    scale = max(abs(val - center) for val in x)
    if scale == 0:
        raise ValueError('Os valores de x devem variar.')

    x_norm = [(val - center) / scale for val in x]
    output = {}

    def calc_sse(real_y, pred_y):
        return sum((r - p) ** 2 for r, p in zip(real_y, pred_y))

    if verbose:
        print('\n--- AJUSTE DE CURVAS (MÍNIMOS QUADRADOS) ---')
        print(f'1. Normalização de x: centro = {center:.6f}, escala = {scale:.6f}')
        print(f'   x normalizado (t) = {_format_vector(x_norm, 4)}')

    # 1. Ajustes polinomiais: Reta (grau 1) e Parábola (grau 2)
    for step_idx, (degree, name) in enumerate([(1, 'reta'), (2, 'parabola')], start=2):
        design = [[val ** j for j in range(degree + 1)] for val in x_norm]
        c = least_squares(design, y)
        pred = [sum(c[j] * (val ** j) for j in range(degree + 1)) for val in x_norm]

        if degree == 1:
            coefficients = [
                c[0] - c[1] * center / scale,
                c[1] / scale,
            ]
        else:
            coefficients = [
                c[0] - c[1] * center / scale + c[2] * (center / scale) ** 2,
                c[1] / scale - 2 * c[2] * center / (scale ** 2),
                c[2] / (scale ** 2),
            ]

        sse = calc_sse(y, pred)
        output[name] = {
            'coefficients': coefficients,
            'predictions': pred,
            'sse': sse,
        }

        if verbose:
            print(f'\n{step_idx}. Ajuste Polinomial — {name.upper()} (grau {degree}):')
            print(f'   Coeficientes normalizados (em t): {_format_vector(c, 6)}')
            print(f'   Coeficientes na escala original (em x): {_format_vector(coefficients, 6)}')
            print(f'   Valores previstos G(x): {_format_vector(pred, 4)}')
            print(f'   Soma dos Erros Quadráticos (EQ): {sse:.8f}')

    # 2. Ajuste exponencial linearizado: ln(y) = alpha + beta * x_norm
    log_y = [math.log(val) for val in y]
    alpha, beta = least_squares([[1, val] for val in x_norm], log_y)

    a_exp = math.exp(alpha - beta * center / scale)
    b_exp = beta / scale
    pred_exp = [math.exp(alpha + beta * val) for val in x_norm]
    sse_exp = calc_sse(y, pred_exp)
    log_sse_exp = calc_sse(log_y, [math.log(val) for val in pred_exp])

    output['exponencial'] = {
        'coefficients': [a_exp, b_exp],
        'predictions': pred_exp,
        'sse': sse_exp,
        'log_sse': log_sse_exp,
    }

    if verbose:
        print('\n4. Ajuste Exponencial Linearizado (ln(y) = alpha + beta*t):')
        print(f'   ln(y) = {_format_vector(log_y, 4)}')
        print(f'   Reta linearizada: alpha = {alpha:.6f}, beta = {beta:.6f}')
        print(f'   Convertendo para y = a * exp(b*x): a = {a_exp:.6f}, b = {b_exp:.6f}')
        print(f'   EQ na escala original = {sse_exp:.8f} (EQ no log = {log_sse_exp:.8f})')

    # 3. Ajuste exponencial direto na escala original
    def profile(exponent):
        shift = max(exponent * val for val in x_norm)
        basis = [math.exp(exponent * val - shift) for val in x_norm]
        amplitude = sum(y_val * b_val for y_val, b_val in zip(y, basis)) / sum(b_val * b_val for b_val in basis)
        predictions = [amplitude * b_val for b_val in basis]
        return calc_sse(y, predictions), predictions, amplitude, shift

    middle, step = beta, 1.0
    for _ in range(60):
        left, right = middle - step, middle + step
        err_left = profile(left)[0]
        err_mid = profile(middle)[0]
        err_right = profile(right)[0]

        if err_mid <= min(err_left, err_right):
            break

        middle = left if err_left < err_right else right
        step *= 2
    else:
        raise ValueError('Não foi possível delimitar o mínimo exponencial.')

    bracket_left, bracket_right = left, right

    golden_ratio = (math.sqrt(5) - 1) / 2
    for _ in range(160):
        u = right - golden_ratio * (right - left)
        v = left + golden_ratio * (right - left)
        if profile(u)[0] < profile(v)[0]:
            right = v
        else:
            left = u

    best_exponent = (left + right) / 2
    error, predictions, amplitude, shift = profile(best_exponent)

    b_direct = best_exponent / scale
    a_direct = math.exp(math.log(amplitude) - shift - b_direct * center)

    output['exponencial_original'] = {
        'coefficients': [a_direct, b_direct],
        'predictions': predictions,
        'sse': error,
    }

    if verbose:
        print('\n5. Ajuste Exponencial Direto (Otimização não linear sem logaritmo):')
        print(f'   Intervalo delimitado para o expoente normalizado: [{bracket_left:.4f}, {bracket_right:.4f}]')
        print(f'   Expoente ótimo (Seção Áurea): {best_exponent:.8f}')
        print(f'   Parâmetros originais: a = {a_direct:.6f}, b = {b_direct:.6f}')
        print(f'   EQ mínimo alcançado = {error:.8f}')

    return output


# --- 4. Integração Numérica ---

def integrate(x, y, verbose=True):
    """Calcula a área sob a curva pelas regras do Trapézio e de Simpson 1/3."""
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('Informe ao menos três pontos para comparar os métodos.')

    finite(list(x) + list(y))

    n = len(x) - 1
    h = (x[-1] - x[0]) / n

    if h <= 0:
        raise ValueError('Os pontos devem ser crescentes e igualmente espaçados.')

    for i in range(n):
        step = x[i + 1] - x[i]
        if step <= 0 or not math.isclose(step, h, rel_tol=1e-9, abs_tol=1e-12 * abs(h)):
            raise ValueError('Os pontos devem ser crescentes e igualmente espaçados.')

    if n % 2 != 0:
        raise ValueError('Simpson 1/3 repetida exige número par de subintervalos.')

    endpoints_trap = (y[0] + y[-1]) / 2
    internal_sum = sum(y[1:-1])
    trapezoid_area = h * (endpoints_trap + internal_sum)

    endpoints_simp = y[0] + y[-1]
    odd_points = y[1:-1:2]
    even_points = y[2:-1:2]
    odd_sum = sum(odd_points)
    even_sum = sum(even_points)
    simpson_area = (h / 3) * (endpoints_simp + 4 * odd_sum + 2 * even_sum)

    if verbose:
        print('\n--- INTEGRAÇÃO NUMÉRICA ---')
        print(f'1. Parâmetros da malha: n = {n} subintervalos | passo h = ({x[-1]} - {x[0]}) / {n} = {h:g}')
        print(f'   Pontos y = {_format_vector(y, 4)}')

        print('\n2. Regra do Trapézio Repetida: I = h * [ (y0 + yn)/2 + soma(internos) ]')
        print(f'   Média dos extremos: ({y[0]} + {y[-1]}) / 2 = {endpoints_trap:.6f}')
        print(f'   Soma dos internos (y1..y{n - 1}): {internal_sum:.6f}')
        print(f'   I_Trapézio = {h:g} * ({endpoints_trap:.6f} + {internal_sum:.6f}) = {trapezoid_area:.8f}')

        print('\n3. Regra de Simpson 1/3 Repetida: I = (h/3) * [ (y0 + yn) + 4*soma(ímpares) + 2*soma(pares) ]')
        print(f'   Soma dos extremos (y0 + y{n}): {endpoints_simp:.6f}')
        print(f'   Pontos ímpares {_format_vector(odd_points, 4)} -> soma = {odd_sum:.6f} (x4 = {4 * odd_sum:.6f})')
        print(f'   Pontos pares   {_format_vector(even_points, 4)} -> soma = {even_sum:.6f} (x2 = {2 * even_sum:.6f})')
        print(f'   I_Simpson = ({h:g} / 3) * ({endpoints_simp:.6f} + {4 * odd_sum:.6f} + {2 * even_sum:.6f}) = {simpson_area:.8f}')

    return {
        'h': h,
        'intervals': n,
        'trapezoid': trapezoid_area,
        'simpson': simpson_area,
    }