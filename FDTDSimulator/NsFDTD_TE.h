#ifndef _NS_FDTD_TE_
#define _NS_FDTD_TE_
#include "FDTD_TE.h"

class NsFDTD_TE: public FDTD_TE{
	typedef FDTD_TE super;
private:
	double R_P, R_M;
public:
	NsFDTD_TE();
	~NsFDTD_TE();
	bool calc();
	void field();
private: 
	complex<double> Dx2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i+1,j+1, t)] + p[index(i+1, j-1, t)] - p[index(i,j+1, t)] - p[index(i, j-1, t)])/2.0;
	};

	complex<double> Dy2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i+1,j+1, t)] + p[index(i-1, j+1, t)] - p[index(i+1,j, t)] - p[index(i-1, j, t)])/2.0;
	};

	void CalcE(){
	//todo i+1 - i �ł͂Ȃ� i - (i-1) j�����l H�̕���+0.5�i��ł��邩��
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
	//�d�E�̌v�ZEx
	for(int i=1; i<mField->getNx()-1; i++)
		for(int j=1; j<mField->getNy()-1; j++)
			EX(i,j, +1) = CEX(i,j)*EX(i,j, 0)
								+ CEXLY(i,j)*( R_P*(HZ(i,j+1, 0) - HZ(i,j, 0))
													+ R_M*Dy2_n(Hz,i,j, 0)
													);
    #ifdef _OPENMP
    #pragma omp for
    #endif
	//�d�E�̌v�ZEy
	for(int i=1; i<mField->getNx()-1; i++)
		for(int j=1; j<mField->getNy()-1; j++)
			EY(i,j, +1) = CEY(i,j)*EY(i,j, 0) 
								- CEYLX(i,j)*( R_P*(HZ(i+1,j, 0) - HZ(i,j, 0))
													+ R_M*Dx2_n(Hz,i,j, 0)
													);
}
	}

	//���E�̌v�Z Hz(i+1/2, j+1/2) -> Hz[i,j]
	//todo i - (i-1)�ł͂Ȃ� i+1 - i ? j�����l H�̕���+0.5�i��ł��邩��
	void CalcH(){
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i=1; i<mField->getNx()-1; i++)
		for(int j=1; j<mField->getNy()-1; j++)
			HZ(i,j, +1) = HZ(i,j, 0)
								-CHZLH(i,j)*( EY(i,j  , +1) - EY(i-1,j, +1) )
								+CHZLH(i,j)*( EX(i  ,j, +1) - EX(i,j-1, +1) );
	}


	void absorbing();
	bool EndTask();

	void ReStart(){
		super::Initialize();
		field();
	}
};
#endif //_NS_FDTD_TE	

			/*
	//todo i+1 - i �ł͂Ȃ� i - (i-1) j�����l H�̕���+0.5�i��ł��邩��
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
	//�d�E�̌v�ZEx
	for(int i=2; i<mField->getNx()-1; i++)
		for(int j=2; j<mField->getNy(); j++)
			EX(i,j, +1) = CEX(i,j)*EX(i,j, 0)
								+ CEXLY(i,j)*( R_P*(HZ(i,j, 0) - HZ(i,j-1, 0))
													+ R_M*Dy2_n(Hz,i,j-1, 0)
													);
    #ifdef _OPENMP
    #pragma omp for
    #endif
	//�d�E�̌v�ZEy
	for(int i=2; i<mField->getNx(); i++)
		for(int j=2; j<mField->getNy()-1; j++)
			EY(i,j, +1) = CEY(i,j)*EY(i,j, 0) 
								- CEYLX(i,j)*( R_P*(HZ(i,j, 0) - HZ(i-1,j, 0))
													+ R_M*Dx2_n(Hz,i-1,j, 0)
													);
}
	}

	//���E�̌v�Z Hz(i+1/2, j+1/2) -> Hz[i,j]
	//todo i - (i-1)�ł͂Ȃ� i+1 - i ? j�����l H�̕���+0.5�i��ł��邩��
	void CalcH(){
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i=1; i<mField->getNx()-1; i++)
		for(int j=1; j<mField->getNy()-1; j++)
			HZ(i,j, +1) = HZ(i,j, 0)
								-CHZLH(i,j)*( EY(i+1,j  , +1) - EY(i,j, +1) )
								+CHZLH(i,j)*( EX(i  ,j+1, +1) - EX(i,j, +1) );
	}
	*/
