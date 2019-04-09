/**
 * @file cgSolver.c
 *
 * @author Gabriel de Souza Barreto - GRR20166812
 *
 * @brief Programa principal do cgSolver
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "utils.h"
#include "likutils.h"


int main(int argc, char *argv[])
{
	LIKWID_MARKER_INIT;

	//! Variáveis:
	FILE *arq;				//! Arquivo de saida

	unsigned int n;			//! Tamanho do sistema
	unsigned int k;			//! Numero de diagonais
	unsigned int simk;		//! Numero de diagonais da matriz simetrica
	unsigned int iter;		//! Numero maximo de iteracoes
	unsigned int numiter;	//!  numero de iteracoes
	double tempo;			//! Variavel para inicio da medicao de tempo
	double ttotal;			//! Tempo gasto no metodo inteiro
	double titer;			//! Tempo medio de cada iteracao do metodo
	double timres;			//! Tempo gasto para calculo do residuo
	double tprep;			//! Tempo para preparo da matriz 
	double w;				//! Parametro pre-condicionador
	double err;				//! Erro aproximado absoluto maximo
	double *md;				//! Vetor de diagonais da matriz de coeficientes A
	double *simmat;			//! Vetor de diagonais da nova matriz simetrica AtA
	double *vb;				//! Vetor b dos termos independentes
	double *atvb;			//! Vetor de novos termos independentes Atb
	double *vetx;			//! Vetor da solução do sistema - [x1..xn]
	double *res;			//! Vetor de resíduo r que é utilizado no método
	double *rresid;			//! Vetor de resíduo do trabalho cuja norma ||r|| deve ser retornada

	char teste;				//! Letra que identifica o precondicionante usado
	char *nome;				//! Nome do arquivo de saida
	//!

	srand(20182);
	if (argc < 10)
	{	
		printf("O programa deve ser executado na forma: \n");
		printf("./cgSolver -n <n> -k <k> -p <ω> -i <i> -e <ε> -o <arquivo_saida> \n");
		return (0);
	}
	n = atoi(argv[2]);
	k = atoi(argv[4]);
	w = atof(argv[6]);
	iter = atoi(argv[8]);
	if ((w >= 1.0) || (w < 0.0))
	{
		printf("O w deve ser 0.0 (sem pre-condicionador) ou 0.0<w<1.0 (pre-condicionador de Jacobi) \n");
		return(0);
	}

	//! Dada uma matriz A k-diagonal, AtA resulta numa matriz (2k-1)-diagonal
	simk = (2 * k) - 1;
	md = malloc(n * k * sizeof(double));
	simmat = malloc(n * simk * sizeof(double));
	vb = malloc((n + 1) * sizeof(double));
	atvb = malloc((n + 1) * sizeof(double));
	vetx = malloc((n + 1) * sizeof(double));
	res = malloc((n + 1) * sizeof(double));
	rresid = malloc((n + 1) * sizeof(double));
	
	//! Obs: Vetor x0 para o inicio das iteracoes preenchido de 0s
	for(unsigned int i = 1; i < (n + 1); i++)
	{
		vetx[i] = 0.0;
	}

	if (argc > 11)
	{
		err = atof(argv[10]);
		nome = argv[12];
		arq = fopen(argv[12], "w");
		#ifdef DEBUG
			printf("O erro lido foi: %-11.7g e deveria ser: %g \n", err, EPS);
			printf("\n");
		#endif
	}
	else
	{
		nome = argv[10];
		arq = fopen(argv[10], "w");
		err = (-201820.0); 
	}
		
	fprintf(arq, "# gsb16 Gabriel Barreto\n");
	fprintf(arq, "#\n");

	criaMatrizA(k, n, md);
	criaVetorB(k, n, vb);
	#ifdef DEBUG
		imprimeAb(k, n, md, vb);
	#endif

	tempo = timestamp();
	preparaMatrizA(k, n, md, simmat);
	preparaVetorB(k, n, md, vb, atvb);
	tprep = timestamp() - tempo;
	#ifdef DEBUG
		imprimeAtAtb(simk, n, simmat, atvb);
	#endif
	
	if (w == 0.0)
	{
		#ifdef DEBUG
			printf("Sera executado Gradiente Conjugado basico! \n");
			printf("\n");
		#endif	
		teste = 'A';
		tempo = timestamp();
		numiter = gradienteConjugado(simk, n, simmat, atvb, vetx, iter, err, res, arq);
		ttotal = timestamp() - tempo;
		titer = (ttotal / numiter);
	}
	else
	{
		#ifdef DEBUG	
			printf("Sera executado Gradiente Conjugado com pre-condicionador de Jacobi! \n");
			printf("\n");
		#endif	
		teste = 'J';
		tempo = timestamp();		
		numiter = gradienteConjugadoPre(simk, n, simmat, atvb, vetx, iter, err, res, arq);
		ttotal = timestamp() - tempo;
		titer = (ttotal / numiter); 
	}		

	#ifdef DEBUG
		imprimeXeSistema(k, n, md, vb, vetx);
	#endif

	tempo = timestamp();
	calculaResiduo(k, n, md, vb, vetx, rresid, arq);
	timres = timestamp() - tempo;

	fprintf(arq, "# Tempo PC: %lf\n", tprep);
	fprintf(arq, "# Tempo iter: %lf\n", titer);
	fprintf(arq, "# Tempo residuo: %lf\n", timres);
	fprintf(arq, "#\n");
	fprintf(arq, "%u\n", n);
	for (unsigned int g = 1; g < (n); ++g)
	{
		fprintf(arq, "%.15g ", vetx[g]);
	}
	fprintf(arq, "%.15g\n", vetx[n]);

	#ifdef DEBUG
		printf("\n");
		printf("O tempo total foi: %-11.7g \n", ttotal);
		printf("O tempo por iteracao foi: %-11.7g \n", titer);
		printf("O tempo para calculo do residuo foi: %-11.7g \n", timres);
		printf("O numero de iteracoes foi: %u \n", numiter);
		printf("\n");
	#endif
	printf("Algoritmo executado com %c (A:SEM, J:JACOBI) e resultados gravados no arquivo %s \n", teste, nome);

	fclose(arq);
	
	free(md);
	free(simmat);
	free(vb);
	free(atvb);
	free(vetx);
	free(res);
	free(rresid);

	LIKWID_MARKER_CLOSE;

	return (0);
}
