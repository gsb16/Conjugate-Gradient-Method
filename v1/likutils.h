/**
 * @file likutils.h
 *
 * @author Gabriel de Souza Barreto - GRR20166812
 *
 * @brief Todos os includes e defines necessários para uso da API do LIKWID
 *
 */

#ifndef _LIKINC
#define  _LIKINC
#include <stdlib.h>
#include <stdio.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif
#endif