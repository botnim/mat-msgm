#include "mex.h"
#include "typeGeneral.h"
#include "MRFEnergy.h"
#include "HotStart.h"

#include "time.h"

#include <string.h>

/*
 * A wrapper for general TRW library by Vladimir Kolmogorov.
 *
 * Usage: 
 *   [x e fmsgs] = Gen_TRWS_wrapper_mex(UTerm, PTerm, method, itr, [[imsgs], [isol]] );
 *
 * Inputs:
 *  UTerm - LxN matrix of unary terms (N-number of variables, L-number of labels)
 *          Each term (col) is a pair [Ei(0), Ei(1), Ei(2), ...].'
 *  PTerm - (2+L^2)xE matrix of pairwise terms (E-number of pair-wise relations (edges))
 *          Each term (col) is of the form [i j Eij(0,0), Eij(0,1),... Eij(1,0), Eij(1,1),...].'
 *          Where i,j are the indices of the variables (first variable is of index 1, last one is N)
 *  method- TRWS ('t'), or BP ('b')
 *  itr   - number of iterations
 *
 *  "Hot Start" options, cannot use both (optional):
 *  imsgs - "hot start messages" (2+L)xE matrix, each column [i j mij(0), mij(1)...]
 *			an initial messages from i to j. (optional)
 *  isol  - "hot start" primal solution 1xN vector of init (1-based) labels per variable 
 *          order of labels must match orders of nodes in UTerm. (optional)
 *
 *
 *
 * Outputs:
 *	x	  - double L vectr of final labels (1-based labels)
 *	e	  -	triplet [energy lowerBound nIter]
 *  fmsgs - final state messages (same format as imsgs). 
 *
 * compiling using:
 * >> ...
	  mex -O -largeArrayDims CXXFLAGS="\$CXXFLAGS -Wno-write-strings"...
	  MRFEnergy.cpp ordering.cpp treeProbabilities.cpp minimize.cpp Gen_TRWS_wrapper_mex.cpp...
	  -output Gen_TRWS_wrapper_mex
 *
 */

inline
void null_fcn(...) {}

#ifdef DEBUG
#define DEBUGmexPrintf mexPrintf
#else
#define DEBUGmexPrintf null_fcn  
#endif

// inputs
enum {
    iU = 0,
    iP = 1,
    iM = 2,
    iT = 3,
	iS = 4,
    nI
};

// outputs
enum {
    oX = 0,
    oE,
	oS, // output state
    nO
};

void my_err_function(char* msg) {
    mexErrMsgTxt(msg);
}


template <typename T, typename REAL>
void Gen_TRWS_wrapper(int nout, mxArray* pout[], int nin, const mxArray* pin[]);


void
mexFunction(
    int nout,
    mxArray* pout[],
    int nin,
    const mxArray* pin[])
{
    if (nin <= iT || nin > nI)
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:nin","Expecting at least %d inputs, and no more than %d", iT, nI);
    
    if (nout == 0)
        return;
    if (nout > nO)
         mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:nout", "Expecting %d outputs", nO);
    
    if ( mxIsComplex(pin[iU]) || mxIsSparse(pin[iU]) )
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:unary_term",
                "Unary term must be full real matrix");
    if ( mxIsComplex(pin[iP]) || mxIsSparse(pin[iP]) )
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:pairwise_term",
                "Pair-wise term must be full real matrix");
    if ( mxIsComplex(pin[iT]) || mxIsSparse(pin[iT]) || mxGetNumberOfElements(pin[iT])!=1)
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:iterations",
                "number of iterations must be full real scalar");
    
    if ( mxGetClassID(pin[iP]) != mxGetClassID(pin[iU]) )
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:energy_terms",
                "Both energy terms must be of the same class");
        
    mwSize L = mxGetM(pin[iU]);
    if ( mxGetNumberOfDimensions(pin[iU]) != 2 ||  mxGetM(pin[iU]) != L )
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:unary_term_size",
                "Unary term must be LxN matrix");
    if ( mxGetNumberOfDimensions(pin[iP]) != 2 || mxGetM(pin[iP]) != (2+L*L) )
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:piarwise_term_size",
                "Unary term must be (2+L^2)xE matrix");
    
    if (mxGetNumberOfElements(pin[iM])!=1 || mxGetClassID(pin[iM])!=mxCHAR_CLASS)
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:method",
                "Method must be a char t/b");
	

    switch (mxGetClassID(pin[iU])) {
        case mxDOUBLE_CLASS:
            // only supported class at the moment
            return Gen_TRWS_wrapper<double, double>(nout, pout, nin, pin);
        case mxINT8_CLASS:
        case mxCHAR_CLASS:
//            return Gen_TRWS_wrapper<char, int>(nout, pout, nin, pin);
        case mxSINGLE_CLASS:
//            return Gen_TRWS_wrapper<float, float>(nout, pout, nin, pin);
        case mxUINT8_CLASS:
//            return Gen_TRWS_wrapper<unsigned char, int>(nout, pout, nin, pin);
        case mxINT16_CLASS:
//            return Gen_TRWS_wrapper<short, int>(nout, pout, nin, pin);
        case mxUINT16_CLASS:
//            return Gen_TRWS_wrapper<unsigned short, int>(nout, pout, nin, pin);
        case mxINT32_CLASS:
//            return Gen_TRWS_wrapper<int, int>(nout, pout, nin, pin);
        case mxUINT32_CLASS:
//            return Gen_TRWS_wrapper<unsigned int, int>(nout, pout, nin, pin);
        case mxINT64_CLASS:
        case mxUINT64_CLASS:
        default:
            mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:energy_class",
                    "Unknown/unsupported class %s",mxGetClassName(pin[iU]));
    } 
    return;
}

template <typename T, typename REAL>
void Gen_TRWS_wrapper(int nout, mxArray* pout[], int nin, const mxArray* pin[])
{ 
    
	MRFEnergy<TypeGeneral>* mrf;
	MRFEnergy<TypeGeneral>::NodeId* nodes;
	MRFEnergy<TypeGeneral>::Options options;
	TypeGeneral::REAL energy, lowerBound;
	unsigned int nIter;

    mwSize N = mxGetN(pin[iU]); // number of nodes/variables
    mwSize L = mxGetM(pin[iU]); // number of labels
    mwSize E = mxGetN(pin[iP]); // number of edges/pairs
    
    T* pU = (T*)mxGetData(pin[iU]);
    T* pP = (T*)mxGetData(pin[iP]);

    // get method
    char method = *(char*)mxGetData(pin[iM]);
    if (method != 't' && method != 'b')
        mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:unknown_method", "Unknown method %c, should be one of t/b", method);

    // get number of iterations
    int iterations = static_cast<int>(mxGetScalar(pin[iT]));
    
	mrf = new MRFEnergy<TypeGeneral>(TypeGeneral::GlobalSize(), my_err_function); // construct with an error message function
	nodes = new MRFEnergy<TypeGeneral>::NodeId[N];

	// construct energy
    // add unary terms
    for ( mwSize ii=0; ii < N ; ii++ ) {
        DEBUGmexPrintf("Adding node %d out of %d\n", ii, N);
        nodes[ii] = mrf->AddNode(TypeGeneral::LocalSize(L), 
                TypeGeneral::NodeData(pU + (L*ii)));
    }

    DEBUGmexPrintf("\n\nDone adding nodes \n\n");
    mwSize step = L*L + 2;
    // add pairwise terms
    for ( mwSize ii=0; ii < E ; ii++ ) {
        DEBUGmexPrintf("Adding edge %d out of %d\n", ii, E);
        
        mrf->AddEdge( nodes[ static_cast<mwSize>(pP[step*ii] - 1) ], 
                nodes[ static_cast<mwSize>(pP[step*ii+1] - 1) ],
                TypeGeneral::EdgeData( TypeGeneral::GENERAL, pP+(step*ii+2) ) );
    }

    DEBUGmexPrintf("\n\nDone adding edges \n\n");
    
	// Hot start options
	if ( nin > iS ) {
		// we have initial state
		if ( mxGetClassID(pin[iS]) != mxGetClassID(pin[iU]) ||
			mxIsSparse(pin[iS]) || 
			mxGetNumberOfDimensions(pin[iS]) != 2 ) {
				mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:hotStart","wrong type for hot start");
		}

		// "hot start messages"
		if ( mxGetM(pin[iS]) == (2+L) ) {

				 mwSize nm = mxGetN(pin[iS]); // number of messages
				 // do hot start!
				 MRFEnergy<TypeGeneral>::NodeId* mIds = new MRFEnergy<TypeGeneral>::NodeId[2*nm];
				 double* messages = new double[L*nm];
				 double *phs = mxGetPr(pin[iS]);

				 for ( unsigned int ei=0; ei < nm ; ei++ ) {
					 mIds[ei*2]     = nodes[ static_cast<mwIndex>(phs[ei*(2+L)])-1 ]; // convert from Matlab's 1-index to 0-index
					 mIds[ei*2 + 1] = nodes[ static_cast<mwIndex>(phs[ei*(2+L)+1])-1 ];

					 memcpy(messages+ei*L, phs + (2+L)*ei + 2, L*sizeof(double));
				 }

				 DEBUGmexPrintf("got hot start messages\n");

				 mrf->HotStart(mIds, messages, nm);

				 DEBUGmexPrintf("ready to go\n");

				 delete[] mIds;
				 delete[] messages;

		} else if ( mxGetNumberOfElements(pin[iS]) == N ) {
			DEBUGmexPrintf("got hot start labels\n");
			
			// copy primal sol to Label array
			MRFEnergy<TypeGeneral>::Label* primal_sol = new MRFEnergy<TypeGeneral>::Label[N];
			
			double* phs = mxGetPr(pin[iS]);

			for ( unsigned int ni(0); ni < N ; ni++ ) {
				primal_sol[ni] = static_cast<MRFEnergy<TypeGeneral>::Label>( phs[ni] - 1 ); // convert 1-based labels to 0-based
			}
			if ( method == 't' ) {
				lowerBound = mrf->PrimalInitTRWS(nodes, primal_sol);
			} else {
				lowerBound = mrf->PrimalInit(nodes, primal_sol);
			}
			
			delete[] primal_sol;

			DEBUGmexPrintf("primal init lower bound = %.4f\n", lowerBound);
		} else {
			mexErrMsgIdAndTxt("Gen_TRWS_wrapper_mex:hotStart","wrong size or type for hot start");
		}
	} // done "hot start"

	options.m_iterMax = iterations; // maximum number of iterations
	options.m_eps = 1e-6;
    if ( method == 't' ) {
        /////////////////////// TRW-S algorithm //////////////////////
        
        nIter = mrf->Minimize_TRW_S(options, lowerBound, energy);
    } else {
        //////////////////////// BP algorithm ////////////////////////

		nIter = mrf->Minimize_BP(options, energy);
    }

	if (nout > oS || nout == 1) {
		DEBUGmexPrintf("Dumping current state\n");

		// dump state
		pout[((nout>oS)?oS:0)] = mxCreateDoubleMatrix(2+L, E , mxREAL);
		double* ps = mxGetPr(pout[((nout>oS)?oS:0)]);

		int* pids = new int[2*E];
		double* pmsgs = new double[L*E]; 
		unsigned int ec = mrf->DumpCurrentState(pids, pmsgs);

		DEBUGmexPrintf("Got ec=%d, E=%d\n", ec, E);

		for ( mwSize ei=0; ei < E ; ei++ ) {
			ps[ei*(2+L)]   = static_cast<double>(pids[2*ei]+1);
			ps[ei*(2+L)+1] = static_cast<double>(pids[2*ei+1]+1);
			memcpy(ps+ei*(2+L)+2, pmsgs+ei*(L), (L)*sizeof(double));
		}
		delete[] pids;   
		delete[] pmsgs;

	}

	if ( nout >= 2 ) {
		// get the energy
		energy = mrf->ComputeSolutionAndEnergy();

		// allocate output
		pout[oX] = mxCreateDoubleMatrix(N, 1, mxREAL);
		double* pX = mxGetPr(pout[oX]);

		// read solution
		for ( mwSize ii=0; ii < N ; ii++ ) {
			pX[ii] = static_cast<double>(mrf->GetSolution(nodes[ii]) + 1); // from 0-based labels to 1-based labels.
		}

		if ( nout > oE ) {
			// energy and bound
			pout[oE] = mxCreateDoubleMatrix(1, 3, mxREAL);
			double* pe = mxGetPr(pout[oE]);

			pe[0] = energy; // energy
			pe[1] = lowerBound; // lower bound
			pe[2] = nIter; // actual number of iterations
		}
	}
	// done
	delete[] nodes;
	delete mrf;
}
    
    
