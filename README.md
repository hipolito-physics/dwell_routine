# Poço duplo finito/Finite double well
Implementação das funções de ondas normalizadas para o problema do poço duplo finito.

Para executar a rotina, basta atentar alguns passos:

### 1 Definir os parâmetros

  Deve-se definir os parâmetros (arquivo config.ini), variando o número quântico ( $n$, onde valores ímpares desta variável representam energias pares e valores pares desta variável representam energias ímpares, respectivamente), a barreira de potencial (variáveis $V_B$ e $V_D$) e a largura do poço (variáveis $a$ e $L$).

  É possível que, para casos limites, pode ocasionar algumas indeterminações, gerando erros no cálculo das energias e consequentemente na visualização das funções de onda.

### 2 Executar a rotina 'eig_energy.py'

  Após definir os parâmetros, deve-se executar a rotina mencionada acima, com intuito de encontrar as energias do poço para a implementação das funções de onda.

### 3 Executar a rotina 'eigen_function.py'

  Executando todos os passos anteriores, é possível visualizar a função de onda com os parâmetros estabelecidos no item 1. Para realizar os diferentes niveis de energia, basta alterar o valor do número quântico, $n$, no arquivo 'config.ini'.
  
### 4 Informações adicionais:
  Caso haja mudança nos parâmetros do poço, deve-se executar novamente a rotina 'eig_energy.py'.
