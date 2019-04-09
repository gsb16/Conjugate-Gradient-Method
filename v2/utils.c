/**
 * @file utils.c
 *
 * @author Gabriel de Souza Barreto - GRR20166812
 *
 * @brief Funções usadas em cgSolver
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "utils.h"
#include "likutils.h"


/**
 * @brief      Função que gera os coeficientes de um sistema linear k-diagonal
 *
 * @param      i   coordenada linha do elemento a ser calculado (0<=i<n)
 * @param      j   coordenada coluna do elemento a ser calculado (0<=j<n)
 * @param      k   quantidade de diagonais da matriz A
 *
 * @return 	   valor gerado para o elemento aij
 */
inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ((i==j)?((double)(k<<1)):(1.0)) * ((double)rand() * invRandMax);
}


/**
 * @brief      Função que gera os termos independentes de um sistema linear k-diagonal
 *
 * @param      k   quantidade de diagonais da matriz A
 *
 * @return 	   valor gerado para o elemento b
 */
inline double generateRandomB( unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ((double)(k<<2)) * ((double)rand() * invRandMax);
}


/**
 * @brief      Função que retorna o índice no array de diagonais de um elemento aij
 *
 * @param      i   coordenada linha do elemento
 * @param      j   coordenada coluna do elemento
 * @param      k   quantidade de diagonais da matriz A
 *
 * @return 	   indice de aij no array de diagonais
 */
inline unsigned int indexMap(unsigned int i, unsigned int j, unsigned int k)
{

	return ((k / 2) + ((i - 1) * k) + j - i);
}


/**
 * @brief	   Função que retorna o valor armazenado no array de diagonais de um termo aij
 *
 * @param      i   coordenada linha do elemento
 * @param      j   coordenada coluna do elemento
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   md  array das diagonais da matriz
 *
 * @return	   valor do termo aij
 */
double obtemValor(unsigned int i, unsigned int j, unsigned int k, unsigned int n, double *md)
{	
	//Teste se o elemento está dentro da matriz foi retirado -> parâmetro n está obsoleto
	if ((j > i ? (j - i) : (i - j)) > (k / 2))
	{
		return (0.0);
	}
	
	return (md[indexMap(i, j, k)]);	
}


/**
 * @brief	   Função que cria a matriz A de coeficientes e armazena as diagonais no array adequado
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   md  array das diagonais da matriz
 */
int criaMatrizA(unsigned int k, unsigned int n, double *md, double *md2)
{
	double termo;

	for (unsigned int i = 1; i <= n; ++i)
	{
	 	for (unsigned int j = 1; j <= n; ++j)
	 	{
	 		if ((j > i ? (j - i) : (i - j)) <= (k / 2))
	 		{
	 			termo = generateRandomA(i, j, k);
	 			md[indexMap(i, j, k)] = termo;
	 			md2[indexMap(j, i, k)] = termo;	 			
	 		}
	 	}
	}
	#ifdef DEBUG
		printf("\n");
	#endif
	return (0);
}


/**
 * @brief	   Função que cria o vetor b de termos independentes do sistema
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   vb  array dos termos independentes b
 */
int criaVetorB(unsigned int k, unsigned int n, double *vb)
{
	double termo;

	for (unsigned int i = 1; i < (n + 1); ++i)
	{
		termo = generateRandomB(k);
		vb[i] = termo;
		#ifdef DEBUG
			// printf("b%u=%-20.15g\n", i, termo);
		#endif
	}
	#ifdef DEBUG
		printf("\n");
	#endif
	return (0);
}


/**
 * @brief	   Função que multiplica AtA criando uma matriz simétrica e positiva definida armazenando em array
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   md  array das diagonais da matriz
 * @param 	   simmat  array das diagonais da matriz AtA resultante
 */
int preparaMatrizA(unsigned int k, unsigned int n, double *md, double *md2, double *simmat)
{
	double a;
	double at;
	double soma;
	unsigned int simk = (2*k) - 1;

	for (unsigned int i = 1; i < (n + 1); ++i)
	{
		for (unsigned int j = 1; j < (n + 1); ++j)
		{
			if ((j > i ? (j - i) : (i - j)) <= (simk / 2))
	 		{
	 			soma = 0.0;
	 			for (unsigned int f = 1; f < (n + 1); ++f)
	 			{
	 				at = obtemValor(i, f, k, n, md2);
	 				a = obtemValor(f, j, k, n, md);
	 				#ifdef DEBUG
	 					// printf("at%u%u=%-20.15g\n", f, i, a);
	 					// printf("a%u%u=%-20.15g\n", f, j, a);
	 				#endif
	 				soma += at * a;
	 			}
	 			simmat[indexMap(i, j, simk)] = soma;
	 			#ifdef DEBUG
	 				// printf("ata%u%u=%-20.15g\n", i, j, soma);
	 			#endif
	 		} 
		}
	}
	#ifdef DEBUG
		printf("\n");
	#endif
	return (0);
}


/**
 * @brief	   Função que multiplica Atb criando um novo vetor de termos independentes 
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   md  array das diagonais da matriz
 * @param 	   vb  array dos termos independentes b
 * @param 	   atvb  array dos novos termos independentes Atb resultantes
 */
int preparaVetorB(unsigned int k, unsigned int n, double *md2, double *vb, double *atvb)
{
	double at;
	double soma;

	for (unsigned int f = 1; f < (n + 1); ++f)
	{
		soma = 0.0;
		for (unsigned int g = 1; g < (n + 1); ++g)
		{
			at = obtemValor(f, g, k, n, md2);
			#ifdef DEBUG
				// printf("at%u%u=%-20.15g\n", g, f, at);
				// printf("b%u=%-20.15g\n", g, vb[g]);
			#endif
			soma += at * vb[g];
		}
		atvb[f] = soma;
		#ifdef DEBUG
			// printf("atb%u=%-20.15g\n", f, soma);
		#endif
	}
	#ifdef DEBUG
		printf("\n");
	#endif
	return (0);
}


/**
 * @brief	   Método Gradiente Conjugado sem pré-condicionante 
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   simmat  array das diagonais da matriz
 * @param 	   atvb  array dos termos independentes b
 * @param 	   vetx  array que contém a solução do sistema
 * @param 	   iter  número máximo de iterações
 * @param 	   err  erro máximo permitido
 * @param 	   res  array para uso no método (r)
 * @param 	   arq  arquivo no qual os resultados serão inseridos
 *
 * @return	   número de iterações que foram feitas
 */
double gradienteConjugado(unsigned int k, unsigned int n, double *simmat, double *atvb,
						 double *vetx, unsigned int iter, double err, double *res, FILE *arq)
{
	//! Variáveis:
	double alpha;			//! s
	double beta;			//! m
	double normx;			//! ||x||
	double relerr;			//! erro relativo atual
	double aux0, aux1;		//! aux, aux1
	double *vetxold;		//! vetor x anterior
	double *vetv;			//! v
	double *vetz;			//! z
	unsigned int numiter;	//! numero de iteracoes
	unsigned int indxmax;	//! indice no qual max(|xatual - xold|)
	double soma;			//! variavel auxiliar nos lacos
	//!

	indxmax = 0;

	//! Algoritmo:
	vetv = malloc((n + 1) * sizeof(double));
	vetz = malloc((n + 1) * sizeof(double));
	vetxold = malloc((n + 1) * sizeof(double));
	memcpy(vetxold, vetx, (n + 1)*sizeof(double));	//! x0 = 0
	memcpy(res, atvb, (n + 1)*sizeof(double));		//! r = b
	memcpy(vetv, atvb, (n + 1)*sizeof(double));		//! v = b

	aux0 = 0.0;
	for(unsigned int g = 1; g < (n + 1); g++) 		//! aux = rtr
	{
		aux0 += res[g] * res[g]; 
	}
	
	numiter = 0;	
	do 									//! para k = 0 : max, faca
	{
		numiter++;
		
		for(unsigned int a = 1; a < (n + 1); a++)	//! z = Av 
		{
			soma = 0.0;
			for(unsigned int b = 1; b < (n + 1); b++)
			{
				soma += (obtemValor(a, b, k, n, simmat)) * vetv[b];
			}
			vetz[a] = soma;
		}
		
		soma = 0.0;
		for(unsigned int g = 1; g < (n + 1); g++) 	//! calcula vtz
		{
			soma += vetv[g] * vetz[g]; 
		}
		//! s = aux/vtz 
		alpha = (fabs(soma) < ZERO) ? 0.0 : (aux0 / soma);

		for(unsigned int g = 1; g < (n + 1); g++) 	//! xk+1 = xk + sv
		{
			vetx[g] = vetxold[g] + (alpha * vetv[g]); 
		}

		for(unsigned int g = 1; g < (n + 1); g++) 	//! r = r - sz
		{
			res[g] = res[g] - (alpha * vetz[g]); 
		}

		aux1 = 0.0;
		for(unsigned int g = 1; g < (n + 1); g++)	//! aux1 = rtr
		{
			aux1 += res[g] * res[g]; 
		}

		normx = 0.0;
		for(unsigned int g = 1; g < (n + 1); g++)	//! calcula ||x||
		{
			if (normx < fabs(vetx[g] - vetxold[g]))
			{
				normx = fabs(vetx[g] - vetxold[g]);
				indxmax = g; 
			}
		}		

		//! m = aux1 / aux
		beta = (fabs(aux0) < ZERO) ? 0.0 : (aux1 / aux0);
		aux0 = aux1;

		for(unsigned int g = 1; g < (n + 1); g++)	//! v = r + mv
		{
			vetv[g] = res[g] + (beta * vetv[g]); 
		}

		#ifdef DEBUG	
			printf("\n");
			printf("||X%02u|| = %-11.7g \n", numiter, normx);
			printf("\n");
		#endif	
		
		fprintf(arq, "# iter %u: %.15g\n", numiter, normx);

		//! xold = x
		memcpy(vetxold, vetx, (n + 1)*sizeof(double));
		//! relerr = max(|Xi - Xi-1|) / Xi
		relerr = (fabs(vetx[indxmax]) < ZERO) ? 0.0 : (normx / fabs(vetx[indxmax]));
		//!
	} while ((numiter < iter) && (relerr > err));
	
	free(vetv);
	free(vetz);
	free(vetxold);

	return (numiter);
}


/**
 * @brief	   Método Gradiente Conjugado com pré-condicionante de Jacobi 
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   simmat  array das diagonais da matriz
 * @param 	   atvb  array dos termos independentes b
 * @param 	   vetx  array que contém a solução do sistema
 * @param 	   iter  número máximo de iterações
 * @param 	   err  erro máximo permitido
 * @param 	   res  array para uso no método (r)
 * @param 	   arq  arquivo no qual os resultados serão inseridos
 *
 * @return	   número de iterações que foram feitas
 */
double gradienteConjugadoPre(unsigned int k, unsigned int n, double *simmat, double *atvb,
						 double *vetx, unsigned int iter, double err, double *res, FILE *arq)
{
	//! Variáveis:
	double alpha;			//! s
	double beta;			//! m
	double normx;			//! ||x||
	double relerr;			//! erro relativo atual
	double aux0, aux1;		//! aux, aux1
	double *vetxold;		//! vetor x anterior
	double *vetv;			//! v
	double *vetz;			//! z
	double *vety;			//! y
	double *zeros;			//! vetor cheio de 0s
	unsigned int numiter;	//! numero de iteracoes
	unsigned int indxmax;	//! indice no qual max(|xatual - xold|)
	double soma;			//! variavel auxiliar nos lacos
	//!

	indxmax = 0;

	//! Algoritmo:
	vetz = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	zeros = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	vetxold = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	memcpy(vetxold, vetx, (n + 1)*sizeof(double));	//! x0 = 0
	memcpy(zeros, vetx, (n + 1)*sizeof(double));	//! preenche vetor zeros de 0s
	memcpy(res, atvb, (n + 1)*sizeof(double));		//! r = b
	vety = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	vetv = aligned_alloc(ALIGNMENT, (((n + 1) * sizeof(double)) + (ALIGNMENT - (((n + 1) * sizeof(double)) % ALIGNMENT))));
	
	aux0 = 0.0;
	for(unsigned int g = 1; g < (n + 1); g++)		//! v = M-1b, y = M-1r 
	{
		vetv[g] = (atvb[g] / obtemValor(g, g, k, n, simmat));
		vety[g] = (res[g] / obtemValor(g, g, k, n, simmat));
		aux0 += vety[g] * res[g];	//! aux = ytr
	}
	//#for(unsigned int g = 1; g < (n + 1); g++)		
	//#{
	//#	aux0 += vety[g] * res[g]; 
	//#}
	
	LIKWID_MARKER_START("op1");
	
	numiter = 0;	
	do 									//! para k = 0 : max, faca
	{
		numiter++;
		
		soma = 0.0;
		//!Unroll and Jam z = Av e soma = vtz
		memcpy(vetz, zeros, (n + 1)*sizeof(double));
		for(unsigned int a = 1; a < (1 + (n + 1) - ((n + 1) % STRIDE)); a+=STRIDE)	//! z = Av
		{
			for(unsigned int b = 1; b < (n + 1); b++)
			{
				for(unsigned int c = 0; c < STRIDE; ++c)
				{
					vetz[a + c] += (obtemValor((a + c), b, k, n, simmat)) * vetv[b];
				}
			}
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				soma += vetz[a + c] * vetv[a + c];	//! calcula vtz
			}
		}

		//!remaider loop
		for(unsigned int a = (1 + (n + 1) - ((n + 1) % STRIDE)); a < (n + 1); ++a)
		{
			for(unsigned int b = 1; b < (n + 1); ++b)
			{
				vetz[a] += (obtemValor(a, b, k, n, simmat)) * vetv[b];
			}
			soma += vetz[a] * vetv[a];
		}

		//#for(unsigned int g = 1; g < (n + 1); g++)	
		//#{
		//#	soma += vetv[g] * vetz[g]; 
		//#}
		//! s = aux/vtz
		alpha = (fabs(soma) < ZERO) ? 0.0 : (aux0 / soma);

		normx = 0.0;
		//!Unroll and Jam do calculo da norma e x = xold
		for(unsigned int g = 1; g < (1 + (n + 1) - ((n + 1) % STRIDE)); g+=STRIDE)	
		{
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				vetx[g + c] += (alpha * vetv[g + c]);			//! xk+1 = xk + sv
				if (normx < fabs(vetx[g + c] - vetxold[g + c]))	//! calcula ||x||
				{
					normx = fabs(vetx[g + c] - vetxold[g + c]);
					indxmax = g + c; 
				}
			} 
		}
		for(unsigned int g = (1 + (n + 1) - ((n + 1) % STRIDE)); g < (n + 1); ++g)	//!remainder loop
		{
			vetx[g] += (alpha * vetv[g]);			//! xk+1 = xk + sv
			if (normx < fabs(vetx[g] - vetxold[g]))	//! calcula ||x||
			{
				normx = fabs(vetx[g] - vetxold[g]);
				indxmax = g; 
			}  
		}

		//#for(unsigned int g = 1; g < (n + 1); g++)	//! r = r - sz
		//#{
		//#	res[g] -= (alpha * vetz[g]);
		//#	vety[g] = (res[g] / obtemValor(g, g, k, n, simmat));	//! y = M-1r 
		//#}
		
		//!Unroll and Jam de r = r - sz e y = M-1r
		for(unsigned int g = 1; g < (1 + (n + 1) - ((n + 1) % STRIDE)); g+=STRIDE)	
		{
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				res[g + c] -= (alpha * vetz[g + c]);
				vety[g + c] = (res[g + c] / obtemValor((g + c), (g + c), k, n, simmat));
			} 
		}
		for(unsigned int g = (1 + (n + 1) - ((n + 1) % STRIDE)); g < (n + 1); ++g)	//!remainder loop
		{
			res[g] -= (alpha * vetz[g]);
			vety[g] = (res[g] / obtemValor(g, g, k, n, simmat));  
		}
		//#for(unsigned int g = 1; g < (n + 1); g++)	
		//#{
		//#	vety[g] = (res[g] / obtemValor(g, g, k, n, simmat)); 
		//#}

		//#for(unsigned int g = 1; g < (n + 1); g++)	
		//#{
		//#	if (normx < fabs(vetx[g] - vetxold[g]))
		//#	{
		//#		normx = fabs(vetx[g] - vetxold[g]);
		//#		indxmax = g; 
		//#	}
		//#}		

		aux1 = 0.0;
		//!Unroll and Jam de aux1 = ytr
		for(unsigned int g = 1; g < (1 + (n + 1) - ((n + 1) % STRIDE)); g+=STRIDE)	
		{
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				aux1 += vety[g + c] * res[g + c];
			} 
		}
		for(unsigned int g = (1 + (n + 1) - ((n + 1) % STRIDE)); g < (n + 1); ++g)	//!remainder loop
		{
			aux1 += vety[g] * res[g];  
		}
		
		//#for(unsigned int g = 1; g < (n + 1); g++)	//! aux1 = ytr
		//#{
		//#	aux1 += vety[g] * res[g]; 
		//#}

		//! m = aux1 / aux
		beta = (fabs(aux0) < ZERO) ? 0.0 : (aux1 / aux0);
		aux0 = aux1;

		//#}for(unsigned int g = 1; g < (n + 1); g++)	//! v = y + mv
		//#}{
		//#}	vetv[g] = vety[g] + (beta * vetv[g]); 
		//#}}

		//!Unroll and Jam de v = y + mv
		for(unsigned int g = 1; g < (1 + (n + 1) - ((n + 1) % STRIDE)); g+=STRIDE)	
		{
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				vetv[g + c] = vety[g + c] + (beta * vetv[g + c]);
			} 
		}
		for(unsigned int g = (1 + (n + 1) - ((n + 1) % STRIDE)); g < (n + 1); ++g)	//!remainder loop
		{
			vetv[g] = vety[g] + (beta * vetv[g]);  
		}

		#ifdef DEBUG
			printf("\n");
			printf("||X%02u|| = %-11.7g \n", numiter, normx);
			printf("\n");
		#endif		

		fprintf(arq, "# iter %u: %.15g\n", numiter, normx);

		//! xold = x
		memcpy(vetxold, vetx, (n + 1)*sizeof(double));
		//! relerr = max(|Xi - Xi-1|) / Xi
		relerr = (fabs(vetx[indxmax]) < ZERO) ? 0.0 : (normx / fabs(vetx[indxmax]));
		//!
	} while ((numiter < iter) && (relerr > err));

	LIKWID_MARKER_STOP("op1");
	
	free(vetv);
	free(vetz);
	free(vety);
	free(vetxold);

	return (numiter);
}


/**
 * @brief	   Função para verificação do tempo 
 *
 * @return	   Tempo
 */
double timestamp(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)((tp.tv_sec*1000.0 + tp.tv_usec/1000.0)/1000.0));
}


/**
 * @brief	   Função que calcula a norma ||r|| do resíduo r = b - Ax 
 *
 * @param      k   quantidade de diagonais da matriz
 * @param 	   n   tamanho da matriz quadrada (nxn)
 * @param 	   md  array das diagonais da matriz
 * @param 	   vb  array dos termos independentes b
 * @param 	   vetx  array que contém a solução do sistema
 * @param 	   rresid  array do resíduo
 * @param 	   arq  arquivo no qual os resultados serão inseridos
 *
 * @return	   número de iterações que foram feitas
 */
double calculaResiduo(unsigned int k, unsigned int n, double *md, double *vb,
					 double *vetx, double *rresid, FILE *arq)
{
	LIKWID_MARKER_START("op2");

	double noreucl;
	double normax;

	memcpy(rresid, vb, (n + 1)*sizeof(double));

	//!Unroll and Jam rresid[]
	for(unsigned int a = 1; a < (1 + (n + 1) - ((n + 1) % STRIDE)); a+=STRIDE)		//! r = b - Ax
	{
		for(unsigned int b = 1; b < (n + 1); ++b)
		{
			for(unsigned int c = 0; c < STRIDE; ++c)
			{
				rresid[a + c] -= (obtemValor((a + c), b, k, n, md)) * vetx[b];
			}
		}
	}

	//!remaider loop
	for(unsigned int a = (1 + (n + 1) - ((n + 1) % STRIDE)); a < (n + 1); ++a)
	{
		for(unsigned int b = 1; b < (n + 1); ++b)
		{
			rresid[a] -= (obtemValor(a, b, k, n, md)) * vetx[b];
		}
	}

	//! calcula ||r|| com norma euclidiana e uso de Unroll and Jam
	noreucl = 0.0;
	for(unsigned int g = 1; g < (1 + (n + 1) - ((n + 1) % STRIDE)); g+=STRIDE)
	{
		for(unsigned int c = 0; c < STRIDE; ++c)
		{
			noreucl += (rresid[g + c] * rresid[g + c]);  
		}
	}
	for(unsigned int g = (1 + (n + 1) - ((n + 1) % STRIDE)); g < (n + 1); ++g)	//!remainder loop
	{
		noreucl += (rresid[g] * rresid[g]);  
	}
	noreucl = sqrt(noreucl);
	
	LIKWID_MARKER_STOP("op2");
	
	//! calcula ||r|| com max
	normax = 0.0;
	for(unsigned int g = 1; g < (n + 1); g++)
	{
		normax = (normax < fabs(rresid[g])) ? fabs(rresid[g]) : normax;  
	}

	#ifdef DEBUG
		printf("\n");
		printf("Valor de residuo: \n");
		printf("Residuo ||r|| com MAX: %+5.5g\n", normax);
		printf("Residuo ||r|| feito com norma euclidiana: %+5.5g\n", noreucl);
		printf("\n");
	#endif

	fprintf(arq, "# residuo: %.15g\n", noreucl);

	return (noreucl);
}


/**
 * @brief	   Função para impressão de A e b em tela (só é chamada com make debug)
 */
void imprimeAb(unsigned int k, unsigned int n, double *md, double *vb)
{
	printf("\n");
	printf("A Matriz gerada foi: \n");
	for (unsigned int i = 1; i < (n + 1); ++i)
	{
		for (unsigned int j = 1; j < (n + 1); ++j)
		{
			printf("%02d%02d=%-11.7g ", i, j, obtemValor(i, j, k, n, md));
		}
		printf("\n");
		printf("\n");
	}
	printf("\n");
	printf("\n");
	printf("O vetor B criado foi: \n");
	printf("[ ");
	for (unsigned int i = 0; i < (n + 1); ++i)
	{
		printf("b%d=%-11.7g ", i, vb[i]);
	}
	printf("]\n");
	printf("\n");
}


/**
 * @brief	   Função para impressão da matriz simétrica e novos b em tela (só é chamada com make debug)
 */
void imprimeAtAtb(unsigned int k, unsigned int n, double *simmat, double *atvb)
{
	printf("\n");
	printf("A Matriz simetrica foi: \n");
	for (unsigned int i = 1; i < (n + 1); ++i)
	{
		for (unsigned int j = 1; j < (n + 1); ++j)
		{
			printf("%02d%02d=%-11.7g ", i, j, obtemValor(i, j, k, n, simmat));
		}
		printf("\n");
		printf("\n");
	}
	printf("\n");
	printf("\n");
	printf("O vetor AtB criado foi: \n");
	printf("[ ");
	for (unsigned int i = 0; i < (n + 1); ++i)
	{
		printf("b%d=%-11.7g ", i, atvb[i]);
	}
	printf("]\n");	
	printf("\n");
}


/**
 * @brief	   Função para impressão do x e "solução em tela (só é chamada com make debug)
 */
void imprimeXeSistema(unsigned int k, unsigned int n, double *md, double *vb, double *vetx)
{	
	double soluc;

	printf("\n");
	printf("O vetor X encontrado foi: \n");
	printf("[ ");
	for (unsigned int i = 0; i < (n + 1); ++i)
	{
		printf("x%d=%-11.7g ", i, vetx[i]);
	}
	printf("]\n");
	printf("\n");
	printf("\n");
	printf("A solucao do sistema ficou: \n");
	for (unsigned int i = 1; i < (n + 1); ++i)
	{
		soluc = 0.0;
		for (unsigned int j = 1; j < (n + 1); ++j)
		{
			soluc += obtemValor(i, j, k, n, md) * vetx[j];
		}
		printf("Linha %u: %-11.7g = %-11.7g \n", i, soluc, vb[i]);
		printf("\n");
	}
	printf("\n");
}
