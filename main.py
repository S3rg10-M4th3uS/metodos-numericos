#!/usr/bin/env python3
"""Atividade de cálculo numérico com menu, sem dependências externas."""
# Importa bibliotecas nativas do Python para ler comandos do terminal, arquivos JSON, matemática e caminhos de pastas
import argparse
import json
import math
from pathlib import Path
# Importa as funções matemáticas criamos no arquivo src/metodos.py
from src.metodos import gaussian, residual, bridge, regression, integrate

# Descobre a pasta onde este arquivo está salvo e carrega os dados dos exercícios do arquivo problemas.json
ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT/'dados/problemas.json').read_text())


# Resolve automaticamente todos os 4 exercícios usando os dados padrão do arquivo JSON
def solve_all(data=DATA):
    # Pega os dados do problema 01.3 (Minas) e resolve o sistema linear por Eliminação de Gauss
    p = data['01.3']
    x = gaussian(p['matrix'], p['rhs'])
    # Retorna um dicionário com os resultados calculados de todos os 4 problemas
    return {'01.3': {'volumes': x, 'residual': residual(p['matrix'], x, p['rhs'])},
            '02.1': bridge(data['02.1']['voltage'], data['02.1']['resistances']),
            '03.3': regression(**data['03.3']), '04.1': integrate(**data['04.1'])}


# Lê uma sequência de números digitados pelo usuário no teclado e valida se estão corretos
def numbers(prompt, count=None):
    # Troca vírgula por ponto, separa os números pelos espaços e converte cada um para decimal (float)
    values = [float(v.replace(',', '.')) for v in input(prompt).split()]
    # Verifica se a lista não está vazia, se tem a quantidade exata esperada e se não há valores infinitos/inválidos
    if not values or (count is not None and len(values) != count) or not all(map(math.isfinite, values)):
        raise ValueError('Quantidade incorreta ou número inválido. Separe os valores por espaços.')
    return values


# Formata e exibe na tela os resultados de cada exercício de forma organizada
def show(key, result):
    print('\nResultado —', key)
    # Se for o problema 01.3, mostra o volume de cada mina, o erro máximo (resíduo) e avisa se algum volume deu negativo
    if key == '01.3':
        for i, value in enumerate(result['volumes'], 1):
            print(f'Mina {i}: {value:.8f} m³')
        print(f"Resíduo máximo: {result['residual']:.3e}")
        if any(v < 0 for v in result['volumes']):
            print('Atenção: volume negativo; solução sem viabilidade física para estas demandas.')
    # Se for o problema 02.1, mostra as tensões VB e VC, as 6 correntes, o total de iterações e o sentido das correntes
    elif key == '02.1':
        print('VB = %.10f V; VC = %.10f V' % tuple(result['voltages']))
        for i, value in enumerate(result['currents'], 1):
            print(f'I{i} = {value:.10f} A')
        print(f"Iterações: {len(result['history'])}; resíduo nodal: {result['residual']:.3e} A")
        print('Sentidos: I1 A→B; I2 A→C; I3 C→D; I4 B→D; I5 B→C; I6 retorno da fonte.')
    # Se for o problema 03.3, monta a fórmula matemática de cada curva ajustada e exibe o Erro Quadrático (EQ)
    elif key == '03.3':
        for name, model in result.items():
            c = model['coefficients']
            formula = f'{c[0]:.10f} exp({c[1]:.10f} x)' if name.startswith('exponencial') else ' + '.join(f'({v:.10f}) x^{i}' for i,v in enumerate(c))
            print(f"{name}: G(x) = {formula}; EQ = {model['sse']:.10f}")
        print('exponencial: ajuste de ln(y); exponencial_original: ajuste direto de y. Todos os EQ na escala original.')
    # Caso contrário (problema 04.1), mostra o tamanho do passo (h), os subintervalos e as áreas pelo Trapézio e por Simpson
    else:
        print(f"h = {result['h']:g}; subintervalos = {result['intervals']}")
        print(f"Trapézio: {result['trapezoid']:.10f} m²; Simpson 1/3: {result['simpson']:.10f} m²")


# Permite que o usuário digite seus próprios dados para testar cada método numérico
def custom(key):
    # Para o problema 01.3, pede o tamanho do sistema, lê linha por linha da matriz A e o vetor b, e resolve por Gauss
    if key == '01.3':
        n = int(input('Número de incógnitas: '))
        if n < 1 or n > 100:
            raise ValueError('Use de 1 a 100 incógnitas.')
        a = [numbers(f'Linha {i+1} de A (coeficientes decimais): ', n) for i in range(n)]
        b = numbers('Vetor b: ', n)
        x = gaussian(a, b)
        return {'volumes': x, 'residual': residual(a,x,b)}
    # Para o problema 02.1, pede a tensão da fonte, os 5 resistores e a tolerância de parada para calcular a ponte
    if key == '02.1':
        e = numbers('Tensão E (V): ', 1)[0]
        r = numbers('R1 R2 R3 R4 R5 (ohms): ', 5)
        tol = numbers('Tolerância absoluta (sugestão 1e-10): ', 1)[0]
        return bridge(e,r,tol)
    # Para os problemas 03.3 e 04.1, lê os valores de x e y (garantindo que profundidades do rio não sejam negativas)
    x = numbers('Valores de x, separados por espaços: ')
    y = numbers('Valores de y/profundidades: ', len(x))
    if key == '04.1' and any(v < 0 for v in y):
        raise ValueError('Profundidades devem ser não negativas.')
    # Chama o ajuste de curvas (03.3) ou o cálculo de área (04.1) conforme a escolha
    return regression(x,y) if key == '03.3' else integrate(x,y)


# Função principal que controla o programa (tanto por comandos diretos no terminal quanto pelo menu interativo)
def main():
    # Configura as opções que podem ser passadas direto ao rodar o script pelo terminal (--todos e --json)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--todos', action='store_true', help='Resolver os quatro enunciados e sair')
    parser.add_argument('--json', type=Path, help='Salvar resultados completos neste arquivo')
    args = parser.parse_args()
    # Se o usuário usou --todos ou --json, resolve tudo direto, mostra na tela, salva o arquivo (se pedido) e encerra
    if args.todos or args.json:
        results = solve_all()
        for key, result in results.items():
            show(key, result)
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(results, ensure_ascii=False, indent=2)+'\n')
        return
    # Mapeia as opções numéricas do menu (1 a 4) para os códigos dos problemas
    keys = {'1':'01.3', '2':'02.1', '3':'03.3', '4':'04.1'}
    # Loop infinito que mantém o menu rodando até o usuário escolher sair
    while True:
        print('\nCÁLCULO NUMÉRICO\n1 — Minas (Gauss)\n2 — Ponte de Wheatstone (Gauss-Seidel)\n3 — Ajustes por mínimos quadrados\n4 — Área do rio\n5 — Resolver todos os enunciados\n0 — Sair')
        try:
            # Lê a opção escolhida e remove espaços extras nas pontas
            choice = input('Opção: ').strip()
            # Se digitar 0, sai do loop e encerra o programa
            if choice == '0':
                break
            # Se digitar 5, resolve e exibe todos os 4 problemas usando os dados do enunciado
            if choice == '5':
                for key, result in solve_all().items():
                    show(key,result)
            # Se digitar de 1 a 4, pergunta se quer usar os dados prontos do enunciado ou digitar novos dados
            elif choice in keys:
                key = keys[choice]
                mode = input('1 — Dados do enunciado; 2 — Novos dados: ').strip()
                if mode not in ('1','2'):
                    raise ValueError('Escolha 1 ou 2.')
                show(key, solve_all()[key] if mode == '1' else custom(key))
            else:
                print('Opção inválida.')
        # Captura erros de digitação ou estouro de cálculo e exibe uma mensagem amigável sem travar o programa
        except (ValueError, OverflowError) as exc:
            print(f'Erro: {exc}')
        # Captura interrupções pelo teclado (como Ctrl+C) para encerrar o programa de forma limpa
        except (EOFError, KeyboardInterrupt):
            print('\nPrograma encerrado.')
            break


# Garante que o menu só seja iniciado se este arquivo for executado diretamente
if __name__ == '__main__':
    main()