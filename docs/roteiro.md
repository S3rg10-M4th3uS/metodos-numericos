# Roteiro de apresentação — aproximadamente 9 minutos

O problema escolhido para exposição é o 04.1. Os outros três estão resolvidos
no relatório e no programa. A divisão abaixo dá a cada integrante cerca de
três minutos, incluindo a demonstração e o encerramento.

## Sérgio — slides 1 a 3 — 3 minutos

Apresentar a equipe e a necessidade de calcular uma seção de rio. Explicar que
as medições fornecem profundidades, e a integral da profundidade em relação à
distância fornece área. Mostrar que a largura é 20 m, que existem 11 pontos
e que, portanto, há 10 subintervalos de 2 m. Destacar profundidade zero nas
margens e a leitura dos rótulos da figura, sem estimar valores pelo desenho.

## Vitor — slide 4 e demonstração do slide 7 — 3 minutos

Explicar a aproximação por segmentos de reta, a soma das áreas dos trapézios
e o peso metade das extremidades. Mostrar a soma interna 31,6 e a área 63,2 m².
Abrir o terminal já posicionado na pasta do projeto. Executar `python3 main.py`,
escolher 4, depois 1. Em seguida, escolher 4, depois 2; fornecer `0 1 2` e
`0 2 0`. Mostrar que o programa calcula novos dados e permite sair com 0.
É possível antecipar o slide 7 e depois retornar ao slide 5.

## Gabrielle — slides 5, 6, 8 e 9 — 3 minutos

Explicar a aproximação por parábolas, os pesos alternados 4 e 2 e a necessidade
de quantidade par de subintervalos. Mostrar as somas 18,2 e 13,4 e o resultado
66,4 m². Comparar com 63,2 m², sem chamar a diferença de erro verdadeiro.
Apresentar brevemente os demais resultados e os testes; concluir com as duas
estimativas. Todos participam das perguntas.

## Perguntas para ensaio

- Por que o espaçamento é 2 m? São 20 m divididos em 10 subintervalos.
- Por que 11 pontos não significa 11 subintervalos? Cada intervalo fica entre dois pontos consecutivos.
- Simpson é sempre melhor? Não; a vantagem de ordem pressupõe suavidade e não determina o erro de um perfil desconhecido.
- Por que Simpson exige número par de subintervalos? Cada parábola integra dois subintervalos.
- Por que a área não é negativa? A profundidade foi definida positiva para baixo; integra-se sua magnitude.
- A diferença de 3,2 m² é o erro? Não conhecemos a área exata para afirmar isso.
- Por que a ponte foi resolvida por tensões? A formulação nodal reduz o sistema e garante diagonal dominante; as correntes seguem da lei de Ohm.
- Por que há duas exponenciais? Minimizar resíduos de logaritmos é diferente de minimizar resíduos dos dados originais.
- Por que a parábola pode extrapolar mal? Seu coeficiente quadrático é positivo; fora da faixa dos dados ela volta a crescer.

## Checklist antes do envio

- Confirmar data da equipe: o PDF diz 17 ou 19/11/2026; mensagem e cabeçalho têm divergências de ano.
- Todos devem estudar os quatro problemas e estar presentes.
- Conferir a divisão atribuída no relatório e adequá-la ao trabalho efetivamente realizado.
- Abrir os PDFs no computador da apresentação e ensaiar a demonstração.
- Publicar o conteúdo de código do ZIP em repositório ou pasta compartilhada; obter e testar o link.
- Enviar relatório e link do código pela plataforma após a apresentação, até 20/11/2026.
- Guardar comprovante de envio. A plataforma e o destino do repositório não foram informados; não houve envio automático.
