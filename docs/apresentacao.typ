#set page(width: 28cm, height: 15.75cm, margin: (x: 1.6cm, y: 0.9cm), fill: rgb("f5f8fc"), footer: context align(right, text(size: 10pt, fill: rgb("516078"))[Cálculo Numérico · #counter(page).display()]))
#set text(font: "DejaVu Sans", size: 17pt, fill: rgb("17283e"), lang: "pt")
#set par(leading: .35em)
#set heading(numbering: none)
#let title(body) = { text(size: 26pt, weight: "bold", fill: rgb("125b83"), body); v(.4cm) }
#let note(body) = block(fill: rgb("e0edf5"), inset: 10pt, radius: 5pt, width: 100%, body)
#title[Da profundidade à área do rio]
#text(size: 25pt)[Problema 04.1 · Integração numérica]
#v(.8cm)
Trapézios e Simpson 1/3 repetidas
#v(1cm)
#text(size: 16pt)[Sérgio Matheus Araujo Silva \
Vitor Massena dos Santos \
Gabrielle Cordeiro Santana]
#v(.4cm)
#text(size: 14pt)[UNIVASF · Cálculo Numérico · Novembro de 2026]
#pagebreak()
#title[1. O que precisamos calcular?]
A área da seção transversal é usada no estudo da vazão e de reservatórios.
#v(.3cm)
Conhecemos profundidades em pontos discretos, mas não uma função exata do fundo.
#v(.4cm)
$ A = integral_0^20 d(x) dif x $
#v(.4cm)
#note[Distância e profundidade em metros → área em metros quadrados.]
#pagebreak()
#title[2. Dados da figura]
#text(size: 15pt)[Fonte: PDF da atividade, página 5. Valores transcritos dos rótulos da figura.]
#v(.4cm)
#table(columns: 6, inset: 9pt, stroke: .5pt + rgb("b3c9db"),
[x (m)], [0], [2], [4], [6], [8],
[d (m)], [0], [1,8], [2], [4], [4],
[x (m)], [10], [12], [14], [16], [18],
[d (m)], [6], [4], [3,6], [3,4], [2,8])
#v(.3cm)
Margem direita: x = 20 m; profundidade = 0 m.
#note[11 pontos · 10 subintervalos · espaçamento h = 2 m]
#pagebreak()
#title[3. Regra dos Trapézios]
Aproximamos o perfil por segmentos de reta entre medições.
#v(.3cm)
$ A_T = h ((d_0 + d_n)/2 + sum_(i=1)^(n-1) d_i) $
#v(.3cm)
Soma das profundidades internas: 31,6 m.
#v(.3cm)
#note[#text(size: 26pt, weight: "bold")[A = 2 × 31,6 = 63,2 m²]]
#pagebreak()
#title[4. Regra de Simpson 1/3]
Aproximamos o perfil por parábolas a cada dois subintervalos.
#v(.3cm)
#text(size: 18pt)[Pesos: 1, 4, 2, 4, 2, 4, 2, 4, 2, 4, 1.]
#v(.4cm)
Soma dos índices ímpares: 18,2 m. \
Soma dos índices pares internos: 13,4 m.
#v(.3cm)
#note[#text(size: 23pt, weight: "bold")[A = (2/3) × (4 × 18,2 + 2 × 13,4) \
A = 66,4 m²]]
#text(size: 15pt)[Requisitos: espaçamento uniforme e número par de subintervalos.]
#pagebreak()
#title[5. Como interpretar os resultados?]
#table(columns: (2fr, 1fr), inset: 10pt, stroke: .5pt + rgb("b3c9db"),
[*Método*], [*Área*], [Trapézios], [63,2 m²], [Simpson 1/3], [66,4 m²])
#v(.3cm)
Diferença: *3,2 m²*, ou *5,06%* do valor dos Trapézios.
#v(.3cm)
#note[A diferença entre métodos não é o erro exato. O perfil contínuo real não foi fornecido.]
#pagebreak()
#title[6. Demonstração do programa]
#text(size: 17pt)[
1. Executar `python3 main.py`.
2. Escolher *4*, depois *1*: reproduzir o enunciado.
3. Escolher *4*, depois *2*: inserir novos dados.
4. Usar x: `0 1 2` e profundidades: `0 2 0`.
5. Conferir: Trapézios = 2 m²; Simpson = 2,6666667 m².
6. Escolher *0* para sair.
]
#v(.3cm)
#note[O exemplo novo mostra que métodos diferentes representam perfis diferentes entre as medições.]
#pagebreak()
#title[7. Os demais problemas da equipe]
#text(size: 18pt)[
*01.3 — Minas:* Gauss com pivoteamento parcial. \
Volumes: 2.416,6667; 9.193,3333; 4.690 m³.

*02.1 — Ponte:* Gauss-Seidel nas tensões nodais. \
10 iterações; corrente total de 0,3586956522 A.

*03.3 — Regressão:* menor erro observado com a parábola. \
EQ = 2,08495696; exponencial direta: EQ = 2,79133895.
]
#v(.4cm)
#note[Os quatro problemas têm código reutilizável, resultados e explicação no relatório.]
#pagebreak()
#title[8. Verificação e participação]
#text(size: 18pt)[
*11 testes automatizados aprovados.*

Sérgio: contexto, dados e modelagem da integral. \
Vitor: regra dos Trapézios e demonstração. \
Gabrielle: Simpson, comparação e conclusão.

Todos: revisão dos quatro problemas e respostas às perguntas.
]
#v(.4cm)
#note[Resultado do rio: *63,2 m²* pelos Trapézios e *66,4 m²* por Simpson.]
