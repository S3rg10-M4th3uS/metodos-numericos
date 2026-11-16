"""Verificação de resultados, propriedades matemáticas e entradas inválidas."""
import math
import unittest
from src.metodos import gaussian, seidel, bridge, regression, integrate
from main import solve_all


class NumericalTests(unittest.TestCase):
    def test_gaussian_pivot_and_preservation(self):
        a, b = [[0,2],[1,3]], [4,7]
        self.assertEqual(gaussian(a,b), [1,2])
        self.assertEqual(a, [[0,2],[1,3]])

    def test_singular(self):
        with self.assertRaises(ValueError):
            gaussian([[1,2],[2,4]],[1,2])

    def test_mine_conservation(self):
        volumes = solve_all()['01.3']['volumes']
        for actual, expected in zip(volumes,[7250/3,27580/3,4690]):
            self.assertAlmostEqual(actual,expected,places=8)
        self.assertAlmostEqual(sum(volumes),16300)

    def test_bridge_kirchhoff(self):
        result = bridge(30,[20,120,120,120,120])
        i1,i2,i3,i4,i5,i6 = result['currents']
        for error in [i1-i4-i5,i2+i5-i3,i6-i3-i4,20*i1+120*i4-30,
                      120*i2+120*i3-30,120*i4-120*i3-120*i5]:
            self.assertAlmostEqual(error,0,places=9)
        self.assertAlmostEqual(i6,33/92,places=9)

    def test_balanced_bridge(self):
        result = bridge(12,[100]*5)
        self.assertAlmostEqual(result['currents'][4],0,places=9)
        self.assertAlmostEqual(result['currents'][5],.12,places=9)

    def test_seidel_failure_and_initial(self):
        with self.assertRaises(ValueError):
            seidel([[1,2],[2,1]],[1,1],max_iter=10)
        with self.assertRaises(ValueError):
            seidel([[0,1],[1,2]],[1,2])
        x,_ = seidel([[4,1],[1,3]],[6,7],initial=[1,1])
        self.assertAlmostEqual(x[0],1)
        self.assertAlmostEqual(x[1],2)

    def test_exact_regressions(self):
        x = [0,1,2,3,4]
        r = regression(x,[2+3*v+v*v for v in x])
        self.assertLess(r['parabola']['sse'],1e-20)
        r = regression(x,[3*math.exp(-.2*v) for v in x])
        for name in ['exponencial','exponencial_original']:
            self.assertLess(r[name]['sse'],1e-14)
            self.assertAlmostEqual(r[name]['coefficients'][0],3,places=6)
        # Resíduos da regressão linear devem ser ortogonais às colunas 1 e x.
        r = regression([0,1,2,3],[1,3,2,5])['reta']
        errors = [u-v for u,v in zip([1,3,2,5],r['predictions'])]
        self.assertAlmostEqual(sum(errors),0)
        self.assertAlmostEqual(sum(i*v for i,v in enumerate(errors)),0)

    def test_regression_validation(self):
        for x,y in [([1,1,1],[1,2,3]),([1,2,3],[1,0,3]),([1,2],[2,3])]:
            with self.assertRaises(ValueError):
                regression(x,y)

    def test_integration_exactness(self):
        x = [0,.5,1,1.5,2]
        self.assertAlmostEqual(integrate(x,[v**3 for v in x])['simpson'],4)
        self.assertAlmostEqual(integrate(x,[2*v+1 for v in x])['trapezoid'],6)
        r = solve_all()['04.1']
        self.assertAlmostEqual(r['trapezoid'],63.2)
        self.assertAlmostEqual(r['simpson'],66.4)

    def test_integration_validation(self):
        for x,y in [([0,1,3],[0,1,0]),([0,1,2,3],[0,1,2,0]),([2,1,0],[0,1,0])]:
            with self.assertRaises(ValueError):
                integrate(x,y)

    def test_invalid_numbers(self):
        with self.assertRaises(ValueError):
            gaussian([[float('nan')]],[1])
        with self.assertRaises(ValueError):
            bridge(30,[0,120,120,120,120])


if __name__ == '__main__':
    unittest.main()
