#!/usr/bin/env python3
"""Atividade de cálculo numérico com menu, sem dependências externas."""
import argparse
import json
import math
from pathlib import Path
from src.metodos import gaussian, residual, bridge, regression, integrate

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT/'dados/problemas.json').read_text())


def solve_all(data=DATA):
    p = data['01.3']
    x = gaussian(p['matrix'], p['rhs'])
    return {'01.3': {'volumes': x, 'residual': residual(p['matrix'], x, p['rhs'])},
            '02.1': bridge(data['02.1']['voltage'], data['02.1']['resistances']),
            '03.3': regression(**data['03.3']), '04.1': integrate(**data['04.1'])}


def numbers(prompt, count=None):
    values = [float(v.replace(',', '.')) for v in input(prompt).split()]
    if not values or (count is not None and len(values) != count) or not all(map(math.isfinite, values)):
        raise ValueError('Quantidade incorreta ou número inválido. Separe os valores por espaços.')
    return values


def show(key, result):
    print('\nResultado —', key)
    if key == '01.3':
        for i, value in enumerate(result['volumes'], 1):
            print(f'Mina {i}: {value:.8f} m³')
        print(f"Resíduo máximo: {result['residual']:.3e}")
        if any(v < 0 for v in result['volumes']):
            print('Atenção: volume negativo; solução sem viabilidade física para estas demandas.')
    elif key == '02.1':
        print('VB = %.10f V; VC = %.10f V' % tuple(result['voltages']))
        for i, value in enumerate(result['currents'], 1):
            print(f'I{i} = {value:.10f} A')
        print(f"Iterações: {len(result['history'])}; resíduo nodal: {result['residual']:.3e} A")
        print('Sentidos: I1 A→B; I2 A→C; I3 C→D; I4 B→D; I5 B→C; I6 retorno da fonte.')
    elif key == '03.3':
        for name, model in result.items():
            c = model['coefficients']
            formula = f'{c[0]:.10f} exp({c[1]:.10f} x)' if name.startswith('exponencial') else ' + '.join(f'({v:.10f}) x^{i}' for i,v in enumerate(c))
            print(f"{name}: G(x) = {formula}; EQ = {model['sse']:.10f}")
        print('exponencial: ajuste de ln(y); exponencial_original: ajuste direto de y. Todos os EQ na escala original.')
    else:
        print(f"h = {result['h']:g}; subintervalos = {result['intervals']}")
        print(f"Trapézio: {result['trapezoid']:.10f} m²; Simpson 1/3: {result['simpson']:.10f} m²")


def custom(key):
    if key == '01.3':
        n = int(input('Número de incógnitas: '))
        if n < 1 or n > 100:
            raise ValueError('Use de 1 a 100 incógnitas.')
        a = [numbers(f'Linha {i+1} de A (coeficientes decimais): ', n) for i in range(n)]
        b = numbers('Vetor b: ', n)
        x = gaussian(a, b)
        return {'volumes': x, 'residual': residual(a,x,b)}
    if key == '02.1':
        e = numbers('Tensão E (V): ', 1)[0]
        r = numbers('R1 R2 R3 R4 R5 (ohms): ', 5)
        tol = numbers('Tolerância absoluta (sugestão 1e-10): ', 1)[0]
        return bridge(e,r,tol)
    x = numbers('Valores de x, separados por espaços: ')
    y = numbers('Valores de y/profundidades: ', len(x))
    if key == '04.1' and any(v < 0 for v in y):
        raise ValueError('Profundidades devem ser não negativas.')
    return regression(x,y) if key == '03.3' else integrate(x,y)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--todos', action='store_true', help='Resolver os quatro enunciados e sair')
    parser.add_argument('--json', type=Path, help='Salvar resultados completos neste arquivo')
    args = parser.parse_args()
    if args.todos or args.json:
        results = solve_all()
        for key, result in results.items():
            show(key, result)
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(results, ensure_ascii=False, indent=2)+'\n')
        return
    keys = {'1':'01.3', '2':'02.1', '3':'03.3', '4':'04.1'}
    while True:
        print('\nCÁLCULO NUMÉRICO\n1 — Minas (Gauss)\n2 — Ponte de Wheatstone (Gauss-Seidel)\n3 — Ajustes por mínimos quadrados\n4 — Área do rio\n5 — Resolver todos os enunciados\n0 — Sair')
        try:
            choice = input('Opção: ').strip()
            if choice == '0':
                break
            if choice == '5':
                for key, result in solve_all().items():
                    show(key,result)
            elif choice in keys:
                key = keys[choice]
                mode = input('1 — Dados do enunciado; 2 — Novos dados: ').strip()
                if mode not in ('1','2'):
                    raise ValueError('Escolha 1 ou 2.')
                show(key, solve_all()[key] if mode == '1' else custom(key))
            else:
                print('Opção inválida.')
        except (ValueError, OverflowError) as exc:
            print(f'Erro: {exc}')
        except (EOFError, KeyboardInterrupt):
            print('\nPrograma encerrado.')
            break


if __name__ == '__main__':
    main()
