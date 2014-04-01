#ifndef _ST_FDTD_TE_
#define _ST_FDTD_TE_
#include "FDTD_TE.h"

class StFDTD_TE: public FDTD_TE{
	typedef FDTD_TE super;
public:
	StFDTD_TE();
	~StFDTD_TE();
	bool calc();
	void field();
private:
	void absorbing();
	void cycle();

	void CalcE(){
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
		//�d�E�̌v�ZEx
		for(int i=1; i<mField->getNpx()-1; i++)
			for(int j=1; j<mField->getNpy()-1; j++)
				EX(i,j) = CEX(i,j)*EX(i,j) 
						+ CEXLY(i,j)*( HZX(i,j+1) - HZX(i,j) + HZY(i,j+1) - HZY(i,j) );
    #ifdef _OPENMP
    #pragma omp for
    #endif
		//�d�E�̌v�ZEy
		for(int i=1; i<mField->getNpx()-1; i++)
			for(int j=1; j<mField->getNpy()-1; j++)
				EY(i,j) = CEY(i,j)*EY(i,j)
					    - CEYLX(i,j)*( HZX(i+1,j)-HZX(i,j) + HZY(i+1,j)-HZY(i,j) );
}
	}

	//�ŊO�ǂ͊��S���̂ł͂��邪���S���C���̂ł͂��邽��, H�͍ŊO�ǂ̌v�Z���K�v? todo
	void CalcH(){
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
		for(int i=1; i<mField->getNpx()-1; i++)
			for(int j=1; j<mField->getNpy()-1; j++)
				HZX(i,j) = CHZX(i,j)*HZX(i,j) - CHZXLX(i,j)*(EY(i,j)-EY(i-1,j) );
    #ifdef _OPENMP
    #pragma omp for
    #endif
		for(int i=1; i<mField->getNpx()-1; i++)
			for(int j=1; j<mField->getNpy()-1; j++)
				HZY(i,j) = CHZY(i,j)*HZY(i,j) + CHZYLY(i,j)*(EX(i,j)-EX(i,j-1) );
}

	}
};
#endif //_ST_FDTD_TE_