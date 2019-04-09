## Trabalho de ci164 2018/2

## Autor
Gabriel de Souza Barreto
GRR20166812

##Arquivos
###Lista de arquivos:
- *cgSolver.c* : Programa principal
- *utils.c*: Funções auxiliares usadas em cgSolver

## Makefile
###Opções:
- *make*: Compila o programa executável cgSolver
- *make clean*: Deleta o executável e a documentação gerada
- *make doc*: Gera a documentação do projeto e link simbólico doc.html
- *make debug*: Compila cgSolver com -DDEBUG e -g para uso do GDB e 
 impressão detalhada de todos os valores na tela durante a execução.

## Execução
O programa deve ser executado na forma:
##### ./cgSolver -n n -k k -p ω -i i -e ε -o arquivo_saida
Onde:
- n: Dimensão do Sistema Linear
- k: Número de diagonais da matriz A (k ímpar)
- ω: Pré-condicionador a ser utilizado (ω=0.0 -> SEM e 0.0<ω<1.0 -> Jacobi) \
- i: Número máximo de iterações a serem executadas
- ε: [OPCIONAL] Erro aproximado absoluto máximo, sendo: 
##### norma máxima (relativa) em x = (max(|xi - xi-1| * 1/|xi|) < ε)
- arquivo_saida: Caminho completo para o arquivo que vai conter a solução

### Biblioteca usada para medição de performance:
 - [LIKWID](https://github.com/RRZE-HPC/likwid)
 
