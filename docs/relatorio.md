---
title: "Cálculo Numérico — Atividade Complementar"
subtitle: "Problemas 01.3, 02.1, 03.3 e 04.1"
author: "Sérgio Matheus Araujo Silva; Vitor Massena dos Santos; Gabrielle Cordeiro Santana"
date: "Novembro de 2026"
lang: pt-BR
---

# Introdução

Este relatório apresenta a modelagem, os métodos, a implementação e a verificação de quatro problemas da atividade complementar de Cálculo Numérico da Fundação Universidade Federal do Vale do São Francisco (UNIVASF): determinação de volumes de mineração, correntes de uma ponte de Wheatstone, regressão por mínimos quadrados e integração da seção transversal de um rio.

Os dados foram transcritos do arquivo `Atividade_complementar2026.2_LES710C.pdf`, incluindo as figuras das páginas 2 e 5. Os algoritmos foram implementados em Python, usando apenas a biblioteca padrão. O programa permite executar os exemplos, fornecer novos dados, realizar novos cálculos e sair pelo menu. Os números exibidos são arredondados apenas na apresentação; os cálculos usam ponto flutuante de dupla precisão.

# Problema 01.3 — Volumes de mineração

## Dados e modelagem

Sejam $x_1$, $x_2$ e $x_3$ os volumes, em m³, extraídos das minas 1, 2 e 3. As composições são:

| Mina | Areia | Cascalho fino | Cascalho grosso |
|:--|--:|--:|--:|
| 1 | 55% | 30% | 15% |
| 2 | 25% | 45% | 30% |
| 3 | 25% | 20% | 55% |

As demandas são 4.800 m³ de areia, 5.800 m³ de cascalho fino e 5.700 m³ de cascalho grosso. Consideram-se as frações volumétricas fornecidas, sem perdas de extração ou processamento. Cada equação representa o balanço de um material:

$$\begin{aligned}
0{,}55x_1+0{,}25x_2+0{,}25x_3&=4800,\\
0{,}30x_1+0{,}45x_2+0{,}20x_3&=5800,\\
0{,}15x_1+0{,}30x_2+0{,}55x_3&=5700.
\end{aligned}$$

As linhas da matriz correspondem aos materiais e as colunas, às minas. As porcentagens são convertidas em frações antes do cálculo.

## Método direto

Aplica-se eliminação de Gauss com pivoteamento parcial. Em cada coluna, escolhe-se o elemento de maior módulo disponível como pivô, trocam-se as linhas quando necessário e eliminam-se os elementos abaixo da diagonal. Para a linha $i$, usa-se o multiplicador $m_{ik}=a_{ik}/a_{kk}$ e a operação $L_i\leftarrow L_i-m_{ik}L_k$. O sistema triangular resultante é resolvido por substituição regressiva:

$$x_i=\frac{b_i-\sum_{j=i+1}^{n}a_{ij}x_j}{a_{ii}}.$$

A implementação rejeita pivôs nulos ou numericamente muito pequenos em relação à escala da matriz. O pivoteamento reduz problemas de arredondamento, mas não elimina o efeito de um sistema mal condicionado.

## Resultado e conferência

| Mina | Volume extraído (m³) |
|:--|--:|
| 1 | 2.416,66666667 |
| 2 | 9.193,33333333 |
| 3 | 4.690,00000000 |
| **Total** | **16.300,00000000** |

A substituição fornece:

| Material | Mina 1 | Mina 2 | Mina 3 | Total (m³) |
|:--|--:|--:|--:|--:|
| Areia | 1.329,166667 | 2.298,333333 | 1.172,500000 | 4.800 |
| Cascalho fino | 725,000000 | 4.137,000000 | 938,000000 | 5.800 |
| Cascalho grosso | 362,500000 | 2.758,000000 | 2.579,500000 | 5.700 |

O resíduo máximo $\|Ax-b\|_\infty$ foi zero na execução registrada em ponto flutuante. Isso não significa ausência universal de erro de arredondamento. A soma dos volumes coincide com a demanda total de 16.300 m³, e todos os volumes são positivos.

# Problema 02.1 — Ponte de Wheatstone

## Convenções e equações de Kirchhoff

A figura da página 2 contém os nós A, B, C e D. Adota-se $V_D=0$ V e $V_A=E=30$ V. Os resistores estão nos ramos AB ($R_1$), AC ($R_2$), CD ($R_3$), BD ($R_4$) e BC ($R_5$). Os sentidos seguem as setas da figura: $I_1$: A→B; $I_2$: A→C; $I_3$: C→D; $I_4$: B→D; $I_5$: B→C; $I_6$: retorno de D à fonte.

Os valores são $R_1=20\ \Omega$ e $R_2=R_3=R_4=R_5=120\ \Omega$. Pela lei dos nós:

$$I_1=I_4+I_5,\qquad I_2+I_5=I_3,\qquad I_6=I_1+I_2=I_3+I_4.$$

Pela lei das malhas:

$$R_1I_1+R_4I_4=E,\quad R_2I_2+R_3I_3=E,\quad R_4I_4-R_3I_3-R_5I_5=0.$$

Para aplicar Gauss-Seidel com convergência garantida neste circuito, resolvem-se as duas tensões nodais desconhecidas. Usando a lei de Ohm nas equações dos nós:

$$\begin{aligned}
\left(\frac1{R_1}+\frac1{R_4}+\frac1{R_5}\right)V_B-\frac{V_C}{R_5}&=\frac E{R_1},\\
-\frac{V_B}{R_5}+\left(\frac1{R_2}+\frac1{R_3}+\frac1{R_5}\right)V_C&=\frac E{R_2}.
\end{aligned}$$

Substituindo os valores e multiplicando as equações por 120:

$$\begin{bmatrix}8&-1\\-1&3\end{bmatrix}\begin{bmatrix}V_B\\V_C\end{bmatrix}=\begin{bmatrix}180\\30\end{bmatrix}.$$

A formulação nodal é equivalente às equações das correntes. O programa trabalha com a matriz de condutâncias, sem multiplicação por 120, de modo que o resíduo corresponde ao desequilíbrio de corrente nos nós, em ampères.

## Método iterativo e convergência

Com $V_B^{(0)}=V_C^{(0)}=0$, as atualizações são:

$$V_B^{(k+1)}=\frac{180+V_C^{(k)}}8,\qquad V_C^{(k+1)}=\frac{30+V_B^{(k+1)}}3.$$

O segundo cálculo usa o valor recém-atualizado de $V_B$, característica do método de Gauss-Seidel. A matriz é estritamente diagonal dominante: $8>1$ e $3>1$. O fator de propagação do erro de $V_C$ é $1/24$, portanto a iteração converge. Para resistências estritamente positivas, a matriz nodal continua estritamente diagonal dominante.

O programa exige simultaneamente variação máxima das tensões menor ou igual a $10^{-10}$ V e resíduo nodal máximo menor ou igual a $10^{-10}$ A. Há limite de 10.000 iterações e detecção de valores não finitos. O enunciado 02.1 não fixa uma tolerância; esse critério foi adotado para fornecer resultados estáveis.

| Iteração | $V_B$ (V) | $V_C$ (V) |
|--:|--:|--:|
| 0 | 0 | 0 |
| 1 | 22,50000000 | 17,50000000 |
| 2 | 24,68750000 | 18,22916667 |
| 3 | 24,77864583 | 18,25954861 |
| 10 | 24,78260870 | 18,26086957 |

## Correntes obtidas

Calculam-se $I_1=(E-V_B)/R_1$, $I_2=(E-V_C)/R_2$, $I_3=V_C/R_3$, $I_4=V_B/R_4$ e $I_5=(V_B-V_C)/R_5$.

| Corrente | Valor (A) | Sentido |
|:--|--:|:--|
| $I_1$ | 0,2608695652 | A→B |
| $I_2$ | 0,0978260870 | A→C |
| $I_3$ | 0,1521739130 | C→D |
| $I_4$ | 0,2065217391 | B→D |
| $I_5$ | 0,0543478261 | B→C |
| $I_6$ | 0,3586956522 | Retorno à fonte |

Foram necessárias 10 iterações. O resíduo nodal máximo foi aproximadamente $5{,}51\times10^{-14}$ A. A soma $I_1+I_2$ coincide com $I_3+I_4$. As quedas de tensão em AB–BD e AC–CD somam 30 V. A ponte não está equilibrada, pois $V_B\ne V_C$ e $I_5\ne0$. Como verificação independente, a solução exata do sistema reduzido é $V_B=570/23$ e $V_C=420/23$ V.

# Problema 03.3 — Ajustes por mínimos quadrados

## Dados e critério

| $x_i$ | 0 | 1,5 | 2,6 | 4,2 | 6 | 8,2 | 10 | 11,4 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| $F(x_i)$ | 18 | 13 | 11 | 9 | 6 | 4 | 2 | 1 |

O erro solicitado é $EQ=\sum_{i=1}^{8}[F(x_i)-G(x_i)]^2$. As unidades de $x$ e $F$ não são especificadas no enunciado; o erro tem unidades de $F$ ao quadrado.

Na reta, $G(x)=a_0+a_1x$; na parábola, $G(x)=a_0+a_1x+a_2x^2$. A matriz de projeto contém as colunas $1,x$ ou $1,x,x^2$. A condição de mínimos quadrados é $X^T(Xa-y)=0$. Em vez de formar explicitamente $X^TX$, o programa usa decomposição QR por Gram-Schmidt modificado com reortogonalização e resolve $Ra=Q^Ty$. Os valores de $x$ são centralizados e escalados antes do ajuste; depois, os coeficientes são convertidos à variável original.

## Ajuste exponencial: duas interpretações

Para $G(x)=ae^{bx}$, a transformação $\ln F=\ln a+bx$ permite uma regressão linear, desde que todos os dados de $F$ sejam positivos. Essa técnica minimiza os resíduos **logarítmicos**; ela não minimiza, em geral, o erro quadrático na escala original pedido no enunciado.

Para eliminar essa ambiguidade, apresentam-se os dois ajustes. No ajuste direto, para cada $b$, a amplitude ótima é

$$a(b)=\frac{\sum_i F(x_i)e^{bx_i}}{\sum_i e^{2bx_i}}.$$

Substitui-se essa expressão no erro e realiza-se uma busca unidimensional por seção áurea, após delimitar um intervalo em torno da estimativa logarítmica. A implementação usa expoentes deslocados para evitar estouro numérico. Para novos conjuntos de dados, essa busca é numérica e não constitui prova de mínimo global; para os dados desta atividade, verificou-se também a condição de estacionariedade e a vizinhança da solução.

## Funções e erros

| Modelo | Função ajustada $G(x)$ | $EQ$ na escala original |
|:--|:--|--:|
| Reta | $15{,}61661009-1{,}38799273x$ | 10,84240061 |
| Parábola | $17{,}22496582-2{,}40595107x+0{,}08842997x^2$ | 2,08495696 |
| Exponencial, ajuste direto | $18{,}00741771e^{-0{,}19153019x}$ | 2,79133895 |
| Exponencial, via logaritmos | $20{,}77096430e^{-0{,}23569035x}$ | 13,95932647 |

A exponencial ajustada diretamente é a resposta para mínimos quadrados na escala original. A versão logarítmica é fornecida como comparação e para a interpretação usual de linearização da exponencial. Todos os erros da tabela são calculados com os dados originais e coeficientes sem arredondamento.

Entre os modelos avaliados, a parábola apresentou o menor erro nos pontos observados. Ela possui três parâmetros, enquanto a reta e a exponencial possuem dois. O menor erro de ajuste não garante melhor extrapolação: a parábola volta a crescer fora do intervalo, e a reta pode prever valores negativos. A exponencial permanece positiva e decrescente para os parâmetros encontrados.

# Problema 04.1 — Área da seção transversal do rio

## Leitura da figura

A figura da página 5 mostra largura total de 20 m, 11 pontos igualmente espaçados e profundidade nula nas margens. Assim, há 10 subintervalos e $h=2$ m. Usam-se os valores numéricos indicados na figura, não estimativas feitas pela altura gráfica do desenho.

| $x$ (m) | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $d$ (m) | 0 | 1,8 | 2 | 4 | 4 | 6 | 4 | 3,6 | 3,4 | 2,8 | 0 |

A área procurada é $A=\int_0^{20}d(x)\,dx$.

## Regra dos Trapézios repetida

$$A_T=h\left[\frac{d_0+d_n}{2}+\sum_{i=1}^{n-1}d_i\right].$$

A soma das profundidades internas é 31,6 m. Portanto, $A_T=2\times31{,}6=\mathbf{63{,}2\ m^2}$. Esse resultado corresponde exatamente à integração do perfil formado por segmentos de reta entre as medições.

## Regra de Simpson 1/3 repetida

Simpson 1/3 requer pontos igualmente espaçados e número par de subintervalos, condições satisfeitas por $n=10$:

$$A_S=\frac h3\left[d_0+d_n+4\sum_{i\text{ ímpar}}d_i+2\sum_{\substack{i\text{ par}\\0<i<n}}d_i\right].$$

A soma nos índices ímpares é $1{,}8+4+6+3{,}6+2{,}8=18{,}2$ m. Nos índices pares internos, $2+4+4+3{,}4=13{,}4$ m. Assim:

$$A_S=\frac23[4(18{,}2)+2(13{,}4)]=\mathbf{66{,}4\ m^2}.$$

A diferença entre os métodos é 3,2 m², equivalente a aproximadamente 5,06% da estimativa dos Trapézios. Essa diferença não é um erro exato, pois o perfil contínuo real é desconhecido. Para funções suficientemente suaves, os erros globais são de ordem $h^2$ para Trapézios e $h^4$ para Simpson; a maior ordem não garante melhor aproximação para qualquer perfil real amostrado.

# Implementação, utilização e verificação

O arquivo `main.py` contém o menu e a apresentação dos resultados. `src/metodos.py` contém os algoritmos. Os dados dos enunciados estão em `dados/problemas.json`. O programa funciona com Python 3, sem instalação de bibliotecas adicionais.

```sh
python3 main.py
python3 main.py --todos
python3 main.py --todos --json resultados/resultados.json
python3 -m unittest discover -s tests -v
```

No menu, seleciona-se um problema e, em seguida, os dados do enunciado ou novos dados. Valores de uma mesma linha são separados por espaços; aceita-se vírgula ou ponto decimal, sem separador de milhar. Por exemplo: `0,55 0,25 0,25`. Após cada cálculo, o menu reaparece. A opção 0 encerra a aplicação.

Para mineração, podem ser fornecidos novos sistemas quadrados. Para a ponte, podem ser alteradas a tensão, as cinco resistências e a tolerância. Nos ajustes, podem ser informados novos pares $(x,y)$, desde que haja pelo menos três pontos, posto suficiente e valores positivos de $y$ para a comparação exponencial. Na integração, novas seções devem apresentar profundidades não negativas, espaçamento uniforme e quantidade par de subintervalos, pois o programa compara os dois métodos.

A implementação informa entradas inválidas, sistemas singulares, resistências não positivas, dados incompatíveis e falhas de convergência. Uma solução algébrica com volumes negativos é sinalizada como fisicamente inviável. O arquivo JSON de resultados inclui todo o histórico de Gauss-Seidel e as previsões dos ajustes.

Foram executados 11 testes automatizados, abrangendo pivoteamento, preservação dos dados, singularidade, conservação dos materiais, leis de Kirchhoff, ponte equilibrada, convergência e falha de Gauss-Seidel, regressões conhecidas, ortogonalidade dos resíduos, exatidão de Simpson para um polinômio cúbico, exatidão dos Trapézios para uma reta e rejeição de entradas inválidas. Todos passaram.

# Divisão igualitária das responsabilidades

As responsabilidades estão distribuídas entre modelagem, implementação, revisão e apresentação dos resultados.

**Sérgio Matheus Araujo Silva:** responsável pelo problema 01.3, pela organização do programa e pelas instruções de execução; participa da revisão do problema 02.1. Na apresentação do problema 04.1, explica o contexto, a transcrição dos dados e a modelagem da integral.

**Vitor Massena dos Santos:** responsável pelo problema 02.1 e pela conferência das equações e dos critérios de convergência; participa da revisão do problema 03.3. Na apresentação do problema 04.1, explica e demonstra a regra dos Trapézios.

**Gabrielle Cordeiro Santana:** responsável pelo problema 03.3 e pela organização dos resultados no relatório; participa da revisão do problema 01.3. Na apresentação do problema 04.1, explica a regra de Simpson e compara os resultados.

**Atividades compartilhadas:** desenvolvimento e conferência do problema 04.1, execução dos testes, revisão final do relatório, ensaio da apresentação e demonstração com novos dados. Todos devem compreender os quatro problemas e estar presentes na apresentação.

# Conclusão

Os volumes determinados atendem às três demandas de materiais. A solução iterativa da ponte converge e satisfaz as leis de Kirchhoff. A parábola apresenta o menor erro quadrático nos dados da regressão, e a distinção entre ajustes exponenciais na escala original e logarítmica foi explicitada. Para a seção do rio, obtêm-se 63,2 m² pelos Trapézios e 66,4 m² por Simpson. O programa permite reproduzir os resultados e resolver novos problemas dentro das hipóteses de cada método.

# Referência

UNIVASF. **Atividade complementar de Cálculo Numérico**, `Atividade_complementar2026.2_LES710C.pdf`, problemas 01.3 (p. 2), 02.1 (pp. 2–3), 03.3 (p. 4) e 04.1 (pp. 4–5).
