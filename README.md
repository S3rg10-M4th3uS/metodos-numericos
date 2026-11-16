# Atividade de Cálculo Numérico — UNIVASF

Equipe: Sérgio Matheus Araujo Silva, Vitor Massena dos Santos e Gabrielle Cordeiro Santana.
Problemas: **01.3, 02.1, 03.3 e 04.1**.

## Executar

Requisito: Python 3.10 ou superior, sem bibliotecas externas.

```sh
python3 main.py
python3 main.py --todos
python3 main.py --todos --json resultados/resultados.json
python3 -m unittest discover -s tests -v
```

No Windows, use `py` no lugar de `python3`, por exemplo: `py main.py`.
Não é necessário instalar pacotes adicionais para executar o programa ou os testes.

O menu permite escolher cada problema, usar dados do enunciado ou novos dados,
voltar a calcular e sair com a opção 0. Separe números por espaços e use ponto
ou vírgula decimal, sem separador de milhar. Exemplo: `0,55 0,25 0,25`.
Resistências devem ser positivas. O ajuste conjunto exige y positivo.
A integração compara Trapézios com Simpson 1/3 e exige espaçamento uniforme,
profundidades não negativas e número par de subintervalos.

## Arquivos

- `main.py`: interface em português e execução dos quatro problemas.
- `src/metodos.py`: algoritmos, com comentários em português.
- `dados/problemas.json`: dados transcritos e conferidos no PDF original.
- `resultados/resultados.json`: resultados, previsões e histórico iterativo.
- `tests/test_metodos.py`: testes matemáticos e validação de entradas.
- `docs/relatorio.md`: fonte editável do relatório.
- `docs/apresentacao.typ`: fonte editável dos slides.
- `docs/roteiro.md`: roteiro de apresentação e perguntas para ensaio.
- `entrega/`: relatório em PDF e DOCX, slides em PDF e pacote ZIP.

## Resultados

| Problema | Resultado |
|---|---|
| 01.3 | Volumes 2.416,66666667; 9.193,33333333; 4.690 m³ |
| 02.1 | I1… I6 = 0,2608695652; 0,0978260870; 0,1521739130; 0,2065217391; 0,0543478261; 0,3586956522 A |
| 03.3 | EQ reta 10,84240061; parábola 2,08495696; exponencial direta 2,79133895; exponencial logarítmica 13,95932647 |
| 04.1 | Trapézios: 63,2 m²; Simpson: 66,4 m² |

A exponencial direta minimiza resíduos em y. A exponencial logarítmica minimiza
resíduos em ln(y). Ambas são apresentadas para esclarecer a interpretação.
A busca exponencial é unidimensional e numérica; para dados arbitrários, não há
garantia geral de mínimo global. A formulação nodal da ponte garante convergência
de Gauss-Seidel para resistências positivas.

## Gerar documentos novamente

Esta etapa é opcional: os documentos prontos estão em `entrega/`.
Para regenerá-los, são necessários Pandoc, Typst e a fonte DejaVu Serif
disponíveis no computador. Essas ferramentas não são necessárias para executar
o programa de cálculo numérico.

```sh
pandoc docs/relatorio.md -o entrega/Relatorio_Calculo_Numerico.pdf --pdf-engine=typst -V mainfont='DejaVu Serif' -V fontsize=11pt -V papersize=a4
pandoc docs/relatorio.md -o entrega/Relatorio_Calculo_Numerico.docx
typst compile docs/apresentacao.typ entrega/Apresentacao_Calculo_Numerico.pdf
```
