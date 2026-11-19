#!/usr/bin/env python3
"""Atividade de cálculo numérico com menu, sem dependências externas."""
import argparse
import json
import math
from pathlib import Path
from src.metodos import gaussian, residual, bridge, regression, integrate

# Carrega o arquivo JSON com os dados padrão dos 4 exercícios
ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / 'dados/problemas.json').read_text())


def solve_one(key, data=DATA, verbose=True):
    """Resolve apenas o exercício especificado pela chave."""
    if key == '01.3':
        mines_data = data['01.3']
        volumes = gaussian(mines_data['matrix'], mines_data['rhs'], verbose=verbose)
        return {
            'volumes': volumes,
            'residual': residual(mines_data['matrix'], volumes, mines_data['rhs']),
        }

    if key == '02.1':
        bridge_data = data['02.1']
        return bridge(bridge_data['voltage'], bridge_data['resistances'], verbose=verbose)

    if key == '03.3':
        return regression(**data['03.3'], verbose=verbose)

    if key == '04.1':
        return integrate(**data['04.1'], verbose=verbose)

    raise KeyError(f'Exercício não encontrado: {key}')


def solve_all(data=DATA, verbose=False):
    """Resolve automaticamente os 4 exercícios usando os dados do dicionário informado."""
    return {
        key: solve_one(key, data=data, verbose=verbose)
        for key in ('01.3', '02.1', '03.3', '04.1')
    }


def numbers(prompt, count=None):
    """Lê números separados por espaço no teclado e valida quantidade e valores."""
    error_msg = 'Quantidade incorreta ou número inválido. Separe os valores por espaços.'

    raw_tokens = input(prompt).split()
    if not raw_tokens:
        raise ValueError(error_msg)

    values = [float(token.replace(',', '.')) for token in raw_tokens]

    if count is not None and len(values) != count:
        raise ValueError(error_msg)

    if not all(math.isfinite(val) for val in values):
        raise ValueError(error_msg)

    return values


def show(key, result):
    """Formata e exibe na tela os resultados de cada exercício."""
    print('\nResultado —', key)

    if key == '01.3':
        for i, value in enumerate(result['volumes'], 1):
            print(f'Mina {i}: {value:.8f} m³')
        print(f"Resíduo máximo: {result['residual']:.3e}")

        if any(v < 0 for v in result['volumes']):
            print('Atenção: volume negativo; solução sem viabilidade física para estas demandas.')

    elif key == '02.1':
        vb, vc = result['voltages']
        print(f'VB = {vb:.10f} V; VC = {vc:.10f} V')

        for i, current in enumerate(result['currents'], 1):
            print(f'I{i} = {current:.10f} A')

        iterations = len(result['history'])
        print(f"Iterações: {iterations}; resíduo nodal: {result['residual']:.3e} A")
        print('Sentidos: I1 A→B; I2 A→C; I3 C→D; I4 B→D; I5 B→C; I6 retorno da fonte.')

    elif key == '03.3':
        for name, model in result.items():
            c = model['coefficients']

            if name.startswith('exponencial'):
                formula = f'{c[0]:.10f} exp({c[1]:.10f} x)'
            else:
                terms = [f'({coef:.10f}) x^{i}' for i, coef in enumerate(c)]
                formula = ' + '.join(terms)

            print(f"{name}: G(x) = {formula}; EQ = {model['sse']:.10f}")

        print('exponencial: ajuste de ln(y); exponencial_original: ajuste direto de y. Todos os EQ na escala original.')

    else:
        print(f"h = {result['h']:g}; subintervalos = {result['intervals']}")
        print(f"Trapézio: {result['trapezoid']:.10f} m²; Simpson 1/3: {result['simpson']:.10f} m²")


def custom(key, verbose=True):
    """Pede novos dados ao usuário no teclado e executa o cálculo do exercício escolhido."""
    if key == '01.3':
        n = int(input('Número de incógnitas: '))
        if n < 1 or n > 100:
            raise ValueError('Use de 1 a 100 incógnitas.')

        matrix_a = [numbers(f'Linha {i + 1} de A (coeficientes decimais): ', n) for i in range(n)]
        vector_b = numbers('Vetor b: ', n)
        volumes = gaussian(matrix_a, vector_b, verbose=verbose)

        return {
            'volumes': volumes,
            'residual': residual(matrix_a, volumes, vector_b),
        }

    if key == '02.1':
        voltage = numbers('Tensão E (V): ', 1)[0]
        resistances = numbers('R1 R2 R3 R4 R5 (ohms): ', 5)
        tol = numbers('Tolerância absoluta (sugestão 1e-10): ', 1)[0]
        return bridge(voltage, resistances, tol, verbose=verbose)

    x = numbers('Valores de x, separados por espaços: ')
    y = numbers('Valores de y/profundidades: ', len(x))

    if key == '04.1' and any(depth < 0 for depth in y):
        raise ValueError('Profundidades devem ser não negativas.')

    if key == '03.3':
        return regression(x, y, verbose=verbose)
    return integrate(x, y, verbose=verbose)


def main():
    """Controla a execução via terminal (--todos / --json) ou pelo menu interativo."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--todos', action='store_true', help='Resolver os quatro enunciados e sair')
    parser.add_argument('--json', type=Path, help='Salvar resultados completos neste arquivo')
    args = parser.parse_args()

    # Modo automático: executa direto sem abrir o menu se --todos ou --json forem passados
    if args.todos or args.json:
        results = solve_all(verbose=False)
        for key, result in results.items():
            show(key, result)

        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            formatted_json = json.dumps(results, ensure_ascii=False, indent=2) + '\n'
            args.json.write_text(formatted_json)
        return

    problem_keys = {
        '1': '01.3',
        '2': '02.1',
        '3': '03.3',
        '4': '04.1',
    }

    menu = (
        '\nCÁLCULO NUMÉRICO\n'
        '1 — Minas (Gauss)\n'
        '2 — Ponte de Wheatstone (Gauss-Seidel)\n'
        '3 — Ajustes por mínimos quadrados\n'
        '4 — Área do rio\n'
        '5 — Resolver todos os enunciados\n'
        '0 — Sair'
    )

    while True:
        print(menu)
        try:
            choice = input('Opção: ').strip()

            if choice == '0':
                break

            if choice == '5':
                for key, result in solve_all(verbose=False).items():
                    show(key, result)

            elif choice in problem_keys:
                key = problem_keys[choice]
                mode = input('1 — Dados do enunciado; 2 — Novos dados: ').strip()

                if mode not in ('1', '2'):
                    raise ValueError('Escolha 1 ou 2.')

                if mode == '1':
                    result = solve_one(key, verbose=True)
                else:
                    result = custom(key, verbose=True)

                show(key, result)

            else:
                print('Opção inválida.')

        except (ValueError, OverflowError) as exc:
            print(f'Erro: {exc}')

        except (EOFError, KeyboardInterrupt):
            print('\nPrograma encerrado.')
            break


if __name__ == '__main__':
    main()