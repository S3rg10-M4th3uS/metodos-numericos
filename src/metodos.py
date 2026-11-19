"""Métodos numéricos implementados com a biblioteca padrão do Python."""
import math


# Garante que nenhum valor em uma lista seja um NaN ou Infinito
def finite(values):
    if not all(math.isfinite(v) for v in values):
        raise ValueError('Todos os valores devem ser finitos.')


# Valida um sistema linear Ax=b
def system(a, b):
    n = len(b)
    # Verifica se b não está vazio, se A tem n linhas e se todas as linhas têm n colunas
    if not n or len(a) != n or any(len(row) != n for row in a):
        raise ValueError('A matriz deve ser quadrada e compatível com b.')
    # Achata a matriz A em uma lista única junto com b e verifica se todos os números são válidos
    finite([v for row in a for v in row] + list(b))
    # Retorna a dimensão n do sistema
    return n


# Calcula o resíduo máximo.
def residual(a, x, b):
    # Para cada linha da matriz A, faz o produto escalar com o vetor solução x,
    # subtrai o valor esperado bi, tira o módulo e retorna o maior erro encontrado.
    return max(abs(sum(v * w for v, w in zip(row, x)) - target)
               for row, target in zip(a, b))


def gaussian(a, b):
    """Eliminação de Gauss com pivoteamento parcial, sem alterar as entradas."""
    # Valida o sistema, cria a matriz aumentada convertendo tudo para float
    # e calcula a escala máxima dos coeficientes para detectar matrizes singulares
    n = system(a, b)
    m = [list(map(float, row)) + [float(target)] for row, target in zip(a, b)]
    scale = max(abs(v) for row in a for v in row)
    # Encontra a linha pivot que tem o maior valor absoluto
    # na coluna k; se esse pivot for muito pequeno, o sistema não
    # tem solução única. Troca a linha atual k com a linha pivot
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(m[i][k]))
        if abs(m[pivot][k]) <= 1e-14 * scale:
            raise ValueError('Sistema singular ou numericamente mal condicionado.')
        m[k], m[pivot] = m[pivot], m[k]
        # Escalonamento:
        # zera os elementos abaixo do pivô na coluna k,
        # subtraindo múltiplos da linha k das linhas inferiores i.
        for i in range(k + 1, n):
            factor = m[i][k] / m[k][k]
            m[i][k] = 0.0
            for j in range(k + 1, n + 1):
                m[i][j] -= factor * m[k][j]
    # Substituição retroativa:
    # da última equação até a primeira, vai isolando e descobrindo cada
    # incógnita com base nas variáveis já descobertas
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][j] * x[j] for j in range(i + 1, n))) / m[i][i]
    finite(x)
    return x


# Resolve o sistema linear por tentativas sucessivas (método iterativo de Gauss-Seidel)
def seidel(a, b, tol=1e-10, max_iter=10000, initial=None):
    # Valida o sistema, a tolerância e garante que a diagonal principal não tenha zeros (evita divisão por zero)
    n = system(a, b)
    finite([tol])
    if tol <= 0 or max_iter < 1 or any(a[i][i] == 0 for i in range(n)):
        raise ValueError('Tolerância e limite positivos; diagonal deve ser não nula.')
    # Define o chute inicial (começa tudo com 0.0 se nenhum chute for informado) e valida o tamanho
    x = [0.0] * n if initial is None else list(initial)
    if len(x) != n:
        raise ValueError('Vetor inicial incompatível.')
    finite(x)
    # Lista para guardar o passo a passo de cada tentativa (iteração)
    history = []
    # Repete o cálculo até atingir o limite máximo de tentativas
    for iteration in range(1, max_iter + 1):
        # Guarda os valores da rodada anterior para comparar a evolução
        old = x.copy()
        # Atualiza cada incógnita x[i] usando imediatamente os valores novos já calculados nesta rodada
        for i in range(n):
            x[i] = (b[i] - sum(a[i][j] * x[j] for j in range(n) if i != j)) / a[i][i]
        finite(x)
        # Calcula a maior mudança entre o valor novo e o antigo (delta) e o erro atual (resíduo)
        delta = max(abs(u - v) for u, v in zip(x, old))
        r = residual(a, x, b)
        # Salva o resumo desta rodada no histórico
        history.append({'iteration': iteration, 'x': x.copy(), 'delta': delta, 'residual': r})
        # Se a mudança e o erro já forem menores ou iguais à tolerância, encontrou a resposta
        if delta <= tol and r <= tol:
            return x, history
    # Se estourar o limite de tentativas sem chegar na precisão desejada, avisa o erro
    raise ValueError(f'Gauss-Seidel não convergiu em {max_iter} iterações.')


# Resolve um circuito elétrico de Ponte de Wheatstone (calcula tensões nos nós e correntes)
def bridge(e, resistances, tol=1e-10):
    # Verifica se a tensão e os 5 resistores são números válidos e maiores que zero
    finite([e] + list(resistances))
    if len(resistances) != 5 or any(r <= 0 for r in resistances):
        raise ValueError('Informe cinco resistências estritamente positivas.')
    r1, r2, r3, r4, r5 = resistances
    # Monta a matriz do sistema (A) e o vetor direito (b) usando as Leis de Kirchhoff nos nós B e C
    a = [[1/r1 + 1/r4 + 1/r5, -1/r5], [-1/r5, 1/r2 + 1/r3 + 1/r5]]
    b = [e/r1, e/r2]
    # Descobre as tensões VB e VC resolvendo o sistema por Gauss-Seidel
    (vb, vc), history = seidel(a, b, tol)
    # Calcula as 5 correntes nos resistores (Lei de Ohm: I = V / R) e soma I1 + I2 para achar a corrente total
    currents = [(e-vb)/r1, (e-vc)/r2, vc/r3, vb/r4, (vb-vc)/r5]
    currents.append(currents[0] + currents[1])
    # Devolve um dicionário organizado com todos os dados e resultados do circuito
    return {'matrix': a, 'rhs': b, 'voltages': [vb, vc], 'currents': currents,
            'history': history, 'residual': residual(a, [vb, vc], b)}


# Encontra a melhor curva que se ajusta aos pontos pelo método dos mínimos quadrados
def least_squares(design, y):
    """Decomposição QR por Gram-Schmidt modificado com reortogonalização."""
    # Pega a quantidade de pontos (m) e de coeficientes (n) e verifica se os dados são compatíveis
    m, n = len(y), len(design[0])
    if m < n or len(design) != m or any(len(row) != n for row in design):
        raise ValueError('Dados insuficientes ou matriz incompatível.')
    finite([v for row in design for v in row] + list(y))
    # Prepara as matrizes Q (colunas perpendiculares entre si) e R (triangular superior)
    q, r = [], [[0.0] * n for _ in range(n)]
    # Percorre cada coluna para torná-la perpendicular (ortogonal) às colunas anteriores
    for j in range(n):
        v = [row[j] for row in design]
        original_norm = math.sqrt(sum(t*t for t in v))
        # Repete o processo duas vezes para evitar erros de arredondamento do computador
        for _ in range(2):
            for i in range(j):
                projection = sum(t*u for t, u in zip(q[i], v))
                r[i][j] += projection
                v = [t - projection*u for t, u in zip(v, q[i])]
        # Calcula o tamanho (norma) do vetor resultante; se for quase zero, os dados são redundantes
        r[j][j] = math.sqrt(sum(t*t for t in v))
        if r[j][j] <= 1e-13 * max(original_norm, 1):
            raise ValueError('Dados sem posto suficiente para o ajuste.')
        # Divide o vetor pelo seu tamanho (normaliza) e guarda na matriz Q
        q.append([t/r[j][j] for t in v])
    # Resolve o sistema triangular final usando Eliminação de Gauss para achar os coeficientes
    return gaussian(r, [sum(t*u for t, u in zip(column, y)) for column in q])


# Ajusta 4 tipos de curvas aos pontos (x, y): reta, parábola, exponencial com log e exponencial direta
def regression(x, y):
    # Valida se há pelo menos 3 pares de pontos e se todos os valores de y são maiores que zero
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('Informe ao menos três pares de dados.')
    finite(list(x) + list(y))
    if any(v <= 0 for v in y):
        raise ValueError('O ajuste exponencial logarítmico exige y > 0.')
    # Centraliza e normaliza x (criando a variável t entre -1 e 1) para evitar números gigantes nas contas
    center = sum(x)/len(x)
    scale = max(abs(v-center) for v in x)
    if scale == 0:
        raise ValueError('Os valores de x devem variar.')
    t = [(v-center)/scale for v in x]
    output = {}
    # Ajusta primeiro uma reta (grau 1) e depois uma parábola (grau 2)
    for degree, name in [(1, 'reta'), (2, 'parabola')]:
        c = least_squares([[v**j for j in range(degree+1)] for v in t], y)
        pred = [sum(c[j]*v**j for j in range(degree+1)) for v in t]
        # Desfaz a normalização de t para devolver os coeficientes na escala original de x
        if degree == 1:
            coefficients = [c[0]-c[1]*center/scale, c[1]/scale]
        else:
            coefficients = [c[0]-c[1]*center/scale+c[2]*(center/scale)**2,
                            c[1]/scale-2*c[2]*center/scale**2, c[2]/scale**2]
        # Guarda os coeficientes, os valores previstos pela curva e a soma dos erros ao quadrado (sse)
        output[name] = {'coefficients': coefficients, 'predictions': pred,
                        'sse': sum((u-v)**2 for u, v in zip(y, pred))}
    # Ajusta a curva exponencial aplicando logaritmo natural em y para transformar a curva em uma reta
    alpha, beta = least_squares([[1, v] for v in t], [math.log(v) for v in y])
    # Desfaz o logaritmo e a normalização para obter os parâmetros reais (a e b) e calcula os erros
    a, b = math.exp(alpha-beta*center/scale), beta/scale
    pred = [math.exp(alpha+beta*v) for v in t]
    output['exponencial'] = {'coefficients': [a, b], 'predictions': pred,
                            'sse': sum((u-v)**2 for u, v in zip(y, pred)),
                            'log_sse': sum((math.log(u)-math.log(v))**2 for u, v in zip(y, pred))}

    # Função auxiliar que calcula o erro da exponencial direto na escala original (sem usar logaritmo)
    def profile(exponent):
        # Usa um deslocamento (shift) para evitar que a exponencial estoure o limite numérico do Python
        shift = max(exponent*v for v in t)
        basis = [math.exp(exponent*v-shift) for v in t]
        amplitude = sum(u*v for u,v in zip(y,basis))/sum(v*v for v in basis)
        predictions = [amplitude*v for v in basis]
        return sum((u-v)**2 for u,v in zip(y,predictions)), predictions, amplitude, shift

    # Etapa 1 da busca: começa pelo chute do ajuste logarítmico e dobra o passo até cercar o ponto de menor erro
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
    # Etapa 2 da busca (Seção Áurea): estreita o intervalo encontrado até achar o expoente com o menor erro
    ratio = (math.sqrt(5)-1)/2
    for _ in range(160):
        u, v = right-ratio*(right-left), left+ratio*(right-left)
        if profile(u)[0] < profile(v)[0]:
            right = v
        else:
            left = u
    # Converte o melhor expoente e amplitude de volta para a escala original e salva no resultado
    exponent = (left+right)/2
    error, predictions, amplitude, shift = profile(exponent)
    b = exponent/scale
    a = math.exp(math.log(amplitude)-shift-b*center)
    output['exponencial_original'] = {'coefficients': [a,b], 'predictions': predictions, 'sse': error}
    return output


# Calcula a área sob uma curva (integral numérica) pelas regras do Trapézio e de Simpson 1/3
def integrate(x, y):
    # Verifica se há pelo menos 3 pontos e se todos os valores são números válidos
    if len(x) != len(y) or len(x) < 3:
        raise ValueError('Informe ao menos três pontos para comparar os métodos.')
    finite(list(x) + list(y))
    # Calcula a quantidade de fatias/intervalos (n) e a largura de cada fatia (h)
    n = len(x)-1
    h = (x[-1]-x[0])/n
    # Garante que os pontos de x estejam em ordem crescente e tenham exatamente o mesmo espaçamento
    if h <= 0 or any(x[i+1] <= x[i] or not math.isclose(x[i+1]-x[i], h, rel_tol=1e-9, abs_tol=1e-12*abs(h)) for i in range(n)):
        raise ValueError('Os pontos devem ser crescentes e igualmente espaçados.')
    # Verifica se o número de intervalos é par (condição obrigatória para a regra de Simpson 1/3)
    if n % 2:
        raise ValueError('Simpson 1/3 repetida exige número par de subintervalos.')
    # Retorna o espaçamento h, o número de intervalos e a área calculada pelos dois métodos
    return {'h': h, 'intervals': n,
            'trapezoid': h*(y[0]/2+sum(y[1:-1])+y[-1]/2),
            'simpson': h/3*(y[0]+y[-1]+4*sum(y[1:-1:2])+2*sum(y[2:-1:2]))}