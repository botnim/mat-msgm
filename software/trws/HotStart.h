#ifndef _HOT_START_H_
#define _HOT_START_H_
/*
 * Code for hot-starting TRW-S
 */

#include "MRFEnergy.h"

// Dump messages
template<class T>
unsigned int
MRFEnergy<T>::DumpCurrentState(int* pmIds, double* pmessages)
{
	if (!m_isEnergyConstructionCompleted)
	{
		CompleteGraphConstruction();
	}


	Node* i;
	MRFEdge *e;
	int k;
	
	// get number of labels
	int nl;
	{
		Vector* M = m_nodeFirst->m_firstForward->m_message.GetMessagePtr();
		nl = M->GetArraySize(m_Kglobal, m_nodeFirst->m_K);
	}

	// get backward messages 
	unsigned int ec(0);
	

	for ( i=m_nodeFirst; i; i=i->m_next)
	{	
		for (e=i->m_firstForward; e; e=e->m_nextForward)
		{
			// messages at this stage are "backward" messages
			pmIds[2*ec+1]   = e->m_head->m_ordering;
			pmIds[2*ec] = e->m_tail->m_ordering;

			
				
			Vector* M = e->m_message.GetMessagePtr();
			for (k=0; k<nl ; k++)
			{
				pmessages[ec*(nl) + k] =  static_cast<double>(M->GetArrayValue(m_Kglobal, i->m_K, k));
			}
			ec++;
			if ( ec > m_edgeNum )
				m_errorFn("something not funny about number of edges");
		}
	}
	return ec;
}

//// Dump messages
//template<class T>
//unsigned int
//MRFEnergy<T>::DumpCurrentState(int* pmIds, double* pmessages)
//{
//	if (!m_isEnergyConstructionCompleted)
//	{
//		CompleteGraphConstruction();
//	}
//
//
//	Node* i;
//	MRFEdge *e;
//	int k;
//	
//	// get number of labels
//	int nl;
//	{
//		Vector* M = m_nodeFirst->m_firstForward->m_message.GetMessagePtr();
//		nl = M->GetArraySize(m_Kglobal, m_nodeFirst->m_K);
//	}
//
//	// get messages
//	unsigned int ec(0), nc(0);
//	Vector* Di = (Vector*) (pmessages + m_edgeNum*nl);
//
//	for ( i=m_nodeFirst; i; i=i->m_next)
//	{	
//		Vector* Di = (Vector*) (pmessages + m_edgeNum*nl + nc*nl);
//
//		Di->Copy(m_Kglobal, i->m_K, &i->m_D);
//		pmIds[2*m_edgeNum + 2*nc]     = i->m_ordering;
//		pmIds[2*m_edgeNum + 2*nc + 1] = i->m_ordering;
//
//
//		for (e=i->m_firstForward; e; e=e->m_nextForward)
//		{
//			Di->Add(m_Kglobal, i->m_K, e->m_message.GetMessagePtr());
//
//			// write nodeIds for this edge
//			pmIds[2*ec]   = e->m_tail->m_ordering;
//			pmIds[2*ec+1] = e->m_head->m_ordering;
//
//			Vector* M = e->m_message.GetMessagePtr();
//			
//			if ( M->GetArraySize(m_Kglobal, i->m_K) != nl )
//				m_errorFn("something not funny about nl and m_K");
//
//			for (k=0; k<nl ; k++)
//			{
//				pmessages[ec*nl + k] =  static_cast<double>(M->GetArrayValue(m_Kglobal, i->m_K, k));
//			}
//
//			ec++;			
//			if ( ec > m_edgeNum )
//				m_errorFn("something not funny about number of edges");
//		}
//		for (e=i->m_firstBackward; e; e=e->m_nextBackward)
//		{
//			Di->Add(m_Kglobal, i->m_K, e->m_message.GetMessagePtr());
//		}
//		Di->ComputeAndSubtractMin(m_Kglobal, i->m_K);
//		nc++;
//		if ( nc > m_nodeNum )
//			m_errorFn("something not funny about number of nodes");
//
//	}
//	return ec;
//}

// start from non-empty messages
template<class T>
void 
MRFEnergy<T>::HotStart(const NodeId* mIds, const double* messages, unsigned int nm)
{
	if (!m_isEnergyConstructionCompleted)
	{
		CompleteGraphConstruction();
	}


	Node *i, *j;
	MRFEdge *e(0);

	// update the messages
	for ( unsigned int ei=0 ; ei < nm ; ei++ ) {		
		i = mIds[2*ei] ;
		j = mIds[2*ei+1];

		bool found = false;	
		for (e=i->m_firstForward; e; e=e->m_nextForward) {
			if ( e->m_head == j ) {
				found = true;
				break;
			}
		}
		if ( ! found ) {
//			mexPrintf("backwards\n");
			for (e=i->m_firstBackward; e; e=e->m_nextBackward) {
//				mexPrintf("Edge %d->%d\n", e->m_tail->m_ordering, e->m_head->m_ordering);
				if ( e->m_tail == j ) {
					found = true;
					break;
				}
			}
		}
	
		// found the relevant edge
		if ( found ) {
			Vector* M = e->m_message.GetMessagePtr();
			int nl = M->GetArraySize(m_Kglobal, i->m_K);

			// mexPrintf("Got edge with %d labels\n", nl);

			for (int k=0; k<nl; k++)
				M->SetArrayValue(m_Kglobal, i->m_K, k, static_cast<REAL>(messages[ei*nl + k]) );
// NOT AN ERROR IF MESSAGE WAS NOT FOUND...
//		} else {
//			mexPrintf("Got edge (%d) %d->%d\n", ei, i->m_ordering, j->m_ordering);
//			m_errorFn("could not find edge");
		}
	}
}




// start from primal solution...
template<class T>
typename T::REAL 
MRFEnergy<T>::PrimalInit(const NodeId* mIds, const Label* primal_sol)
{
	if (!m_isEnergyConstructionCompleted)
	{
		CompleteGraphConstruction();
	}

	// make sure all messages are init to zero...
	ZeroMessages();

	Node *i, *j;
	MRFEdge *e(0);
	
	Vector* Di = (Vector*) m_buf;

	// set solution at all nodes
	for ( unsigned int ni=0 ; ni < m_nodeNum ; ni++ )
		SetSol(mIds[ni], primal_sol[ni]);



	////////////////////////////////////////////////
	//               backward pass                //
	////////////////////////////////////////////////

	for (j=m_nodeLast; j; j=j->m_prev)
	{
		// set message from j to i, for i ordered BEFORE j
		// according to GIVEN label of j: lj=l
		// to be M_{ji;k} = \theta_{ji;lk}

		
		// pass messages from i to nodes with smaller m_ordering
		for (e=j->m_firstBackward; e; e=e->m_nextBackward)
		{
			assert(j == e->m_head);
			i = e->m_tail;

			Vector* pV = (Vector*)e->m_message.GetMessagePtr();

			// // If dir==0, then sets dest[kj] += V(ksource,kj).
			// // If dir==1, then sets dest[ki] += V(ki,ksource).
			// // If Swap() has been called odd number of times, then the meaning of dir is swapped.
			// void AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir);
			e->m_message.AddColumn(m_Kglobal, //GlobalSize
				j->m_K, // Ksource
				i->m_K, // Kdest
				j->m_solution, // ksource
				pV, // dest
				0); // dir
			
			
			// subtruct min from message...
			// pV->ComputeAndSubtractMin(m_Kglobal, j->m_K);

		}
	}
	return 0;
}

// start from primal solution...
template<class T>
typename T::REAL 
MRFEnergy<T>::PrimalInitTRWS(const NodeId* mIds, const Label* primal_sol)
{
	if (!m_isEnergyConstructionCompleted)
	{
		CompleteGraphConstruction();
	}

	// make sure all messages are init to zero...
	ZeroMessages();

	// make sure we have weights of messages (tree probabilities)
	SetMonotonicTrees();

	// set solution at all nodes
	for ( unsigned int ni=0 ; ni < m_nodeNum ; ni++ )
		SetSol(mIds[ni], primal_sol[ni]);

	Node *i, *j;
	MRFEdge *e(0);

	Vector* Di = (Vector*) m_buf;
	void* buf = (void*) (m_buf + m_vectorMaxSizeInBytes);

	////////////////////////////////////////////////
	//               backward pass                //
	////////////////////////////////////////////////

	for (i=m_nodeLast; i; i=i->m_prev)
	{
		Di->Copy(m_Kglobal, i->m_K, &i->m_D);
		for (e=i->m_firstBackward; e; e=e->m_nextBackward)
		{
			Di->Add(m_Kglobal, i->m_K, e->m_message.GetMessagePtr());
		}
		for (e=i->m_firstForward; e; e=e->m_nextForward)
		{
			Di->Add(m_Kglobal, i->m_K, e->m_message.GetMessagePtr());
		}

		// set Di[li] = 0
		// Di->SetZero(m_Kglobal, i->m_solution);
		Di->SetArrayValue(m_Kglobal, i->m_K, i->m_solution, 0);

		// set messages from i to nodes with smaller m_ordering
		for (e=i->m_firstForward; e; e=e->m_nextForward)
		{
			assert(e->m_tail == i);
			j = e->m_head;

			Vector* Mji = e->m_message.GetMessagePtr();

			for (int k=0; k<Di->GetArraySize(m_Kglobal, i->m_K); k++) 
			{
				Mji->SetArrayValue(m_Kglobal, i->m_K, k, - e->m_gammaForward * Di->GetArrayValue(m_Kglobal, i->m_K, k));
			}

		}

	}
	return 0;
}

#endif // _HOT_START_H_
