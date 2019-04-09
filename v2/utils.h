/**
 * @file utils.h
 *
 * @author Gabriel de Souza Barreto - GRR20166812
 *
 * @brief Header das funções usadas em cgSolver
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#define ZERO 1.0e-20
#define EPS 1.0e-8
#define ALIGNMENT 16
#define STRIDE 8


double generateRandomA(unsigned int i, unsigned int j, unsigned int k);

double generateRandomB(unsigned int k);

unsigned int indexMap(unsigned int i, unsigned int j, unsigned int k);

double obtemValor(unsigned int i, unsigned int j, unsigned int k, unsigned int n, double *md);

int criaMatrizA(unsigned int k, unsigned int n, double *md, double *md2);

int criaVetorB(unsigned int k, unsigned int n, double *vb);

int preparaMatrizA(unsigned int k, unsigned int n, double *md, double *md2, double *simmat);

int preparaVetorB(unsigned int k, unsigned int n, double *md, double *vb, double *atvb);

double gradienteConjugado(unsigned int k, unsigned int n, double *md, double *atvb,
						 double *vetx, unsigned int iter, double err, double *res, FILE *arq);

double gradienteConjugadoPre(unsigned int k, unsigned int n, double *md, double *atvb,
						 double *vetx, unsigned int iter, double err, double *res, FILE *arq);

double calculaResiduo(unsigned int k, unsigned int n, double *md, double *vb,
					 double *vetx, double *rresid, FILE *arq);

double timestamp(void);

void imprimeAb(unsigned int k, unsigned int n, double *md, double *vb);

void imprimeAtAtb(unsigned int k, unsigned int n, double *simmat, double *atvb);

void imprimeXeSistema(unsigned int k, unsigned int n, double *md, double *vb, double *vetx);
