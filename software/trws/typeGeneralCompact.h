/******************************************************************
typeGeneralCompact.h

Energy function with general interactions:
   E(x)   =   \sum_i D_i(x_i)   +   \sum_ij V_ij(x_i,x_j)
   where x_i \in {0, 1, ..., Ki-1}
   V_ij(ki, kj) are given either as matrices Ki*Kj, or as Potts terms 
   (V_ij(ki, kj) = 0 if ki==kj, and lambda_ij otherwise, with non-negative lambda_ij).

   Inefficient! If possible, use other type*.h files.

*****************************************************************
this implementation is similar to typeGeneral but saves memory in the sense it DOES NOT COPY tables provided
to the energy. 
DO NOT CORRUPT THE DATA for this implementation!!!!
*****************************************************************

organization of matrix V (pair-wise interaction)

 i\j___kj_
  |       |
ki|       | Ki (num labels in i)
  |_______|

  V(ki,kj) = V[ki + Ki*kj] - column stack of V

*****************************************************************

Example usage:

Minimize function E(x,y) = Dx(x) + Dy(y) + lambda*[x != y] + Dz(z) + V(y,z) where 
  x,y \in {0,1,2}, z \in {0,1}
  Dx(0) = 0, Dx(1) = 1, Dx(2) = 2,
  Dy(0) = 3, Dy(1) = 4, Dy(2) = 5,
  lambda = 6,
  [.] is 1 if it's argument is true, and 0 otherwise.
  Dz(0) = 7, Dz(1) = 8,
  V(y,z) = y*y + z



#include <stdio.h>
#include "MRFEnergy.h"

void testGeneral()
{
	MRFEnergy<typeGeneralCompact>* mrf;
	MRFEnergy<typeGeneralCompact>::NodeId* nodes;
	MRFEnergy<typeGeneralCompact>::Options options;
	typeGeneralCompact::REAL energy, lowerBound;

	const int nodeNum = 3; // number of nodes
	typeGeneralCompact::REAL Dx[3]; 
	typeGeneralCompact::REAL Dy[3];
	typeGeneralCompact::REAL Dz[2];
	typeGeneralCompact::REAL V[3*2]; // need different arrays for the tables - do not override !! (BAGON)
	int x, y, z;

	mrf = new MRFEnergy<typeGeneralCompact>(typeGeneralCompact::GlobalSize());
	nodes = new MRFEnergy<typeGeneralCompact>::NodeId[nodeNum];

	// construct energy
	Dx[0] = 0; Dx[1] = 1; Dx[2] = 2;
	nodes[0] = mrf->AddNode(typeGeneralCompact::LocalSize(3), typeGeneralCompact::NodeData(Dx));
	Dy[0] = 3; Dy[1] = 4; Dy[2] = 5;
	nodes[1] = mrf->AddNode(typeGeneralCompact::LocalSize(3), typeGeneralCompact::NodeData(Dy));
	mrf->AddEdge(nodes[0], nodes[1], typeGeneralCompact::EdgeData(typeGeneralCompact::POTTS, 6));
	Dz[0] = 7; Dz[1] = 8;
	nodes[2] = mrf->AddNode(typeGeneralCompact::LocalSize(2), typeGeneralCompact::NodeData(Dz));
	for (y=0; y<3; y++)
	{
		for (z=0; z<2; z++)
		{
			V[y + z*3] = y*y + z;
		}
	}
	mrf->AddEdge(nodes[1], nodes[2], typeGeneralCompact::EdgeData(typeGeneralCompact::GENERAL, V));

	// Function below is optional - it may help if, for example, nodes are added in a random order
	// mrf->SetAutomaticOrdering();

	/////////////////////// TRW-S algorithm //////////////////////
	options.m_iterMax = 30; // maximum number of iterations
	mrf->Minimize_TRW_S(options, lowerBound, energy);

	// read solution
	x = mrf->GetSolution(nodes[0]);
	y = mrf->GetSolution(nodes[1]);

	printf("Solution: %d %d\n", x, y);

	//////////////////////// BP algorithm ////////////////////////
	mrf->ZeroMessages(); // in general not necessary - it may be faster to start 
	                     // with messages computed in previous iterations

	options.m_iterMax = 30; // maximum number of iterations
	mrf->Minimize_BP(options, energy);

	// read solution
	x = mrf->GetSolution(nodes[0]);
	y = mrf->GetSolution(nodes[1]);

	printf("Solution: %d %d\n", x, y);

	// done
	delete nodes;
	delete mrf;
}

*******************************************************************/












#ifndef __typeGeneralCompact_H__
#define __typeGeneralCompact_H__

#include <string.h>
#include <assert.h>


template <class T> class MRFEnergy;


class typeGeneralCompact
{
private:
	struct Vector; // node parameters and messages
	struct Edge; // stores edge information and either forward or backward message

public:
	typedef enum
	{
		GENERAL, // edge information is stored as Ki*Kj matrix. Inefficient!
		POTTS    // edge information is stored as one number (lambdaPotts).
	} Type;

	// types declarations
	typedef int Label;
	typedef double REAL;
	struct GlobalSize; // global information about number of labels
	struct LocalSize; // local information about number of labels (stored at each node)
	struct NodeData; // argument to MRFEnergy::AddNode()
	struct EdgeData; // argument to MRFEnergy::AddEdge()


	struct GlobalSize
	{
	};

	struct LocalSize // number of labels is stored at MRFEnergy::m_Kglobal
	{
		LocalSize(int K);

	private:
	friend struct Vector;
	friend struct Edge;
		int		m_K; // number of labels
	};

	struct NodeData
	{
		NodeData(REAL* data); // data = pointer to array of size MRFEnergy::m_Kglobal

	private:
	friend struct Vector;
	friend struct Edge;
		REAL*		m_data;
	};

	struct EdgeData
	{
		EdgeData(Type type, REAL lambdaPotts); // type must be POTTS
		EdgeData(Type type, REAL* data); // type must be GENERAL. data = pointer to array of size Ki*Kj
		                                 // such that V(ki,kj) = data[ki + Ki*kj]

	private:
	friend struct Vector;
	friend struct Edge;
		Type		m_type;
		union
		{
			REAL	m_lambdaPotts;
			REAL*	m_dataGeneral;
		};
	};







	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////// Visible only to MRFEnergy /////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////

private:
friend class MRFEnergy<typeGeneralCompact>;

	struct Vector
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize K); // returns -1 if invalid K's
		void Initialize(GlobalSize Kglobal, LocalSize K, NodeData data);  // called once when user adds a node
		void Add(GlobalSize Kglobal, LocalSize K, NodeData data); // called once when user calls MRFEnergy::AddNodeData()

		void SetZero(GlobalSize Kglobal, LocalSize K);                            // set this[k] = 0
		void Copy(GlobalSize Kglobal, LocalSize K, Vector* V);                    // set this[k] = V[k]
		void Add(GlobalSize Kglobal, LocalSize K, Vector* V);                     // set this[k] = this[k] + V[k]
		REAL GetValue(GlobalSize Kglobal, LocalSize K, Label k);                  // return this[k]
		REAL ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin);            // return min_k { this[k] }, set kMin
		REAL ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K);              // same as previous, but additionally set this[k] -= vMin (and kMin is not returned)

		static int GetArraySize(GlobalSize Kglobal, LocalSize K);
		REAL GetArrayValue(GlobalSize Kglobal, LocalSize K, int k); // note: k is an integer in [0..GetArraySize()-1].
		                                                            // For Potts functions GetArrayValue() and GetValue() are the same,
		                                                            // but they are different for, say, 2-dimensional labels.
		void SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x);

	private:
	friend struct Edge;
		REAL		m_data[1]; // actual size is MRFEnergy::m_Kglobal
	};

	struct Edge
	{
		static int GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data); // returns -1 if invalid data
		static int GetBufSizeInBytes(int vectorMaxSizeInBytes); // returns size of buffer need for UpdateMessage()
		void Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj); // called once when user adds an edge
		Vector* GetMessagePtr();
		void Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj); // if the client calls this function, then the meaning of 'dir'
								                                               // in distance transform functions is swapped

		// When UpdateMessage() is called, edge contains message from dest to source.
		// The function must replace it with the message from source to dest.
		// The update rule is given below assuming that source corresponds to tail (i) and dest corresponds
		// to head (j) (which is the case if dir==0).
		//
		// 1. Compute Di[ki] = gamma*source[ki] - message[ki].  (Note: message = message from j to i).
		// 2. Compute distance transform: set
		//       message[kj] = min_{ki} (Di[ki] + V(ki,kj)). (Note: message = message from i to j).
		// 3. Compute vMin = min_{kj} m_message[kj].
		// 4. Set m_message[kj] -= vMin.
		// 5. Return vMin.
		//
		// If dir==1 then source corresponds to j, sink corresponds to i. Then the update rule is
		//
		// 1. Compute Dj[kj] = gamma*source[kj] - message[kj].  (Note: message = message from i to j).
		// 2. Compute distance transform: set
		//       message[ki] = min_{kj} (Dj[kj] + V(ki,kj)). (Note: message = message from j to i).
		// 3. Compute vMin = min_{ki} m_message[ki].
		// 4. Set m_message[ki] -= vMin.
		// 5. Return vMin.
		//
		// If Edge::Swap has been called odd number of times, then the meaning of dir is swapped.
		//
		// Vector 'source' must not be modified. Function may use 'buf' as a temporary storage.
		REAL UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* buf);

		// If dir==0, then sets dest[kj] += V(ksource,kj).
		// If dir==1, then sets dest[ki] += V(ki,ksource).
		// If Swap() has been called odd number of times, then the meaning of dir is swapped.
		void AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir);

	protected:

		Type		m_type;

		// message
		Vector*		m_message;
	};

	struct EdgePotts : Edge
	{
	private:
	friend struct Edge;
		REAL	m_lambdaPotts;
	};

	struct EdgeGeneral : Edge
	{
	private:
	friend struct Edge;
		int		m_dir; // 0 if Swap() was called even number of times, 1 otherwise
		REAL*	m_data; // pointer to array of size Ki*Kj STORED EXTERNALLY !!! <<------ (bagon)
	};
};




//////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Implementation ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


inline typeGeneralCompact::LocalSize::LocalSize(int K)
{
	m_K = K;
}

///////////////////// NodeData and EdgeData ///////////////////////

inline typeGeneralCompact::NodeData::NodeData(REAL* data)
{
	m_data = data;
}

inline typeGeneralCompact::EdgeData::EdgeData(Type type, REAL lambdaPotts)
{
	assert(type == POTTS);
	m_type = type;
	m_lambdaPotts = lambdaPotts;
}

inline typeGeneralCompact::EdgeData::EdgeData(Type type, REAL* data)
{
	assert(type == GENERAL);
	m_type = type;
	m_dataGeneral = data;
}

///////////////////// Vector ///////////////////////

inline int typeGeneralCompact::Vector::GetSizeInBytes(GlobalSize Kglobal, LocalSize K)
{
	if (K.m_K < 1)
	{
		return -1;
	}
	return K.m_K*sizeof(REAL);
}
inline void typeGeneralCompact::Vector::Initialize(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	memcpy(m_data, data.m_data, K.m_K*sizeof(REAL));
}

inline void typeGeneralCompact::Vector::Add(GlobalSize Kglobal, LocalSize K, NodeData data)
{
	for (int k=0; k<K.m_K; k++)
	{
		m_data[k] += data.m_data[k];
	}
}

inline void typeGeneralCompact::Vector::SetZero(GlobalSize Kglobal, LocalSize K)
{
	memset(m_data, 0, K.m_K*sizeof(REAL));
}

inline void typeGeneralCompact::Vector::Copy(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	memcpy(m_data, V->m_data, K.m_K*sizeof(REAL));
}

inline void typeGeneralCompact::Vector::Add(GlobalSize Kglobal, LocalSize K, Vector* V)
{
	for (int k=0; k<K.m_K; k++)
	{
		m_data[k] += V->m_data[k];
	}
}

inline typeGeneralCompact::REAL typeGeneralCompact::Vector::GetValue(GlobalSize Kglobal, LocalSize K, Label k)
{
	assert(k>=0 && k<K.m_K);
	return m_data[k];
}

inline typeGeneralCompact::REAL typeGeneralCompact::Vector::ComputeMin(GlobalSize Kglobal, LocalSize K, Label& kMin)
{
	REAL vMin = m_data[0];
	kMin = 0;
	for (int k=1; k<K.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
			kMin = k;
		}
	}

	return vMin;
}

inline typeGeneralCompact::REAL typeGeneralCompact::Vector::ComputeAndSubtractMin(GlobalSize Kglobal, LocalSize K)
{
	REAL vMin = m_data[0];
	for (int k=1; k<K.m_K; k++)
	{
		if (vMin > m_data[k])
		{
			vMin = m_data[k];
		}
	}
	for (int k=0; k<K.m_K; k++)
	{
		m_data[k] -= vMin;
	}

	return vMin;
}

inline int typeGeneralCompact::Vector::GetArraySize(GlobalSize Kglobal, LocalSize K)
{
	return K.m_K;
}

inline typeGeneralCompact::REAL typeGeneralCompact::Vector::GetArrayValue(GlobalSize Kglobal, LocalSize K, int k)
{
	assert(k>=0 && k<K.m_K);
	return m_data[k];
}

inline void typeGeneralCompact::Vector::SetArrayValue(GlobalSize Kglobal, LocalSize K, int k, REAL x)
{
	assert(k>=0 && k<K.m_K);
	m_data[k] = x;
}

///////////////////// EdgeDataAndMessage implementation /////////////////////////

inline int typeGeneralCompact::Edge::GetSizeInBytes(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data)
{
	int messageSizeInBytes = ((Ki.m_K > Kj.m_K) ? Ki.m_K : Kj.m_K)*sizeof(REAL);

	switch (data.m_type)
	{
		case POTTS:
			if (Ki.m_K != Kj.m_K || data.m_lambdaPotts < 0)
			{
				return -1;
			}
			return sizeof(EdgePotts) + messageSizeInBytes;
		case GENERAL:
			return sizeof(EdgeGeneral) + messageSizeInBytes; //BAGON
			// return sizeof(EdgeGeneral) - sizeof(REAL) + Ki.m_K*Kj.m_K*sizeof(REAL) + messageSizeInBytes;
		default:
			return -1;
	}
}

inline int typeGeneralCompact::Edge::GetBufSizeInBytes(int vectorMaxSizeInBytes)
{
	return vectorMaxSizeInBytes;
}

inline void typeGeneralCompact::Edge::Initialize(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj, EdgeData data, Vector* Di, Vector* Dj)
{
	m_type = data.m_type;

	switch (m_type)
	{
		case POTTS:
			((EdgePotts*)this)->m_lambdaPotts = data.m_lambdaPotts;
			m_message = (Vector*)((char*)this + sizeof(EdgePotts));
			break;
		case GENERAL:
			((EdgeGeneral*)this)->m_dir = 0;
			((EdgeGeneral*)this)->m_data = data.m_dataGeneral; // BAGON
			// memcpy(((EdgeGeneral*)this)->m_data, data.m_dataGeneral, Ki.m_K*Kj.m_K*sizeof(REAL));

			m_message = (Vector*)((char*)this + sizeof(EdgeGeneral)); // BAGON
			// m_message = (Vector*)((char*)this + sizeof(EdgeGeneral) - sizeof(REAL) + Ki.m_K*Kj.m_K*sizeof(REAL));
			break;
		default:
			assert(0);
	}

	memset(m_message->m_data, 0, ((Ki.m_K > Kj.m_K) ? Ki.m_K : Kj.m_K)*sizeof(REAL));
}

inline typeGeneralCompact::Vector* typeGeneralCompact::Edge::GetMessagePtr()
{
	return m_message;
}

inline void typeGeneralCompact::Edge::Swap(GlobalSize Kglobal, LocalSize Ki, LocalSize Kj)
{
	if (m_type == GENERAL)
	{
		((EdgeGeneral*)this)->m_dir = 1 - ((EdgeGeneral*)this)->m_dir;
	}
}

inline typeGeneralCompact::REAL typeGeneralCompact::Edge::UpdateMessage(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Vector* source, REAL gamma, int dir, void* _buf)
{
	Vector* buf = (Vector*) _buf;
	REAL vMin;

	if (m_type == POTTS)
	{
		assert(Ksource.m_K == Kdest.m_K);

		int k, kMin;

		m_message->m_data[0] = gamma*source->m_data[0] - m_message->m_data[0];
		kMin = 0;
		vMin = m_message->m_data[0];

		for (k=1; k<Ksource.m_K; k++)
		{
			m_message->m_data[k] = gamma*source->m_data[k] - m_message->m_data[k];
			kMin = 0;
			vMin = buf->m_data[0];
			if (vMin > m_message->m_data[k])
			{
				kMin = k;
				vMin = m_message->m_data[k];
			}
		}

		for (k=0; k<Ksource.m_K; k++)
		{
			m_message->m_data[k] -= vMin;
			if (m_message->m_data[k] > ((EdgePotts*)this)->m_lambdaPotts)
			{
				m_message->m_data[k] = ((EdgePotts*)this)->m_lambdaPotts;
			}
		}
	}
	else if (m_type == GENERAL)
	{
		int ksource, kdest;
		REAL* data = ((EdgeGeneral*)this)->m_data;

		for (ksource=0; ksource<Ksource.m_K; ksource++)
		{
			buf->m_data[ksource] = gamma*source->m_data[ksource] - m_message->m_data[ksource];
		}

		if (dir == ((EdgeGeneral*)this)->m_dir)
		{
			for (kdest=0; kdest<Kdest.m_K; kdest++)
			{
				vMin = buf->m_data[0] + data[0 + kdest*Ksource.m_K];
				for (ksource=1; ksource<Ksource.m_K; ksource++)
				{
					if (vMin > buf->m_data[ksource] + data[ksource + kdest*Ksource.m_K])
					{
						vMin = buf->m_data[ksource] + data[ksource + kdest*Ksource.m_K];
					}
				}
				m_message->m_data[kdest] = vMin;
			}
		}
		else
		{
			for (kdest=0; kdest<Kdest.m_K; kdest++)
			{
				vMin = buf->m_data[0] + data[kdest + 0*Kdest.m_K];
				for (ksource=1; ksource<Ksource.m_K; ksource++)
				{
					if (vMin > buf->m_data[ksource] + data[kdest + ksource*Kdest.m_K])
					{
						vMin = buf->m_data[ksource] + data[kdest + ksource*Kdest.m_K];
					}
				}
				m_message->m_data[kdest] = vMin;
			}
		}

		vMin = m_message->m_data[0];
		for (kdest=1; kdest<Kdest.m_K; kdest++)
		{
			if (vMin > m_message->m_data[kdest])
			{
				vMin = m_message->m_data[kdest];
			}
		}

		for (kdest=0; kdest<Kdest.m_K; kdest++)
		{
			m_message->m_data[kdest] -= vMin;
		}
	}
	else
	{
		assert(0);
	}

	return vMin;
}

inline void typeGeneralCompact::Edge::AddColumn(GlobalSize Kglobal, LocalSize Ksource, LocalSize Kdest, Label ksource, Vector* dest, int dir)
{
	assert(ksource>=0 && ksource<Ksource.m_K);

	int k;

	if (m_type == POTTS)
	{
		for (k=0; k<ksource; k++)
		{
			dest->m_data[k] += ((EdgePotts*)this)->m_lambdaPotts;
		}
		for (k++; k<Kdest.m_K; k++)
		{
			dest->m_data[k] += ((EdgePotts*)this)->m_lambdaPotts;
		}
	}
	else if (m_type == GENERAL)
	{
		REAL* data = ((EdgeGeneral*)this)->m_data;

		if (dir == ((EdgeGeneral*)this)->m_dir)
		{
			for (k=0; k<Kdest.m_K; k++)
			{
				dest->m_data[k] += data[ksource + k*Ksource.m_K];
			}
		}
		else
		{
			for (k=0; k<Kdest.m_K; k++)
			{
				dest->m_data[k] += data[k + ksource*Kdest.m_K];
			}
		}
	}
	else
	{
		assert(0);
	}
}

#endif
