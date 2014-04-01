#ifndef _ST_FDTD_TM_H
#define _ST_FDTD_TM_H
#include "FDTD_TM.h"

class StFDTD_TM: public FDTD_TM{
	typedef FDTD_TM super;
public:
	StFDTD_TM();
	~StFDTD_TM();
	bool calc();
	void field();
private:
	void absorbing();
	void cycle();


	//todo ���E�ߖT��S-FDTD���g���ĂȂ�
	void CalcE(){
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=1; i<mField->getNpx()-1; i++)
			for(int j=1; j<mField->getNpy()-1; j++)
				EZX(i,j) = CEZX(i,j)*EZX(i,j) + CEZXLX(i,j)*(HY(i,j) - HY(i-1,j));

		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=1; i<mField->getNpx()-1; i++)
			for(int j=1; j<mField->getNpy()-1; j++)
				EZY(i,j) = CEZY(i,j)*EZY(i,j) - CEZYLY(i,j)*(HX(i,j) - HX(i,j-1));
		}
	}

	//�z�����E��Ez�ɂ����K�p���Ȃ�����,H�͗̈�̒[�����ʂɌv�Z����(�ł��镪��)
	void CalcH(){	//todo �v�Z�̈� i=0?, 1? j=0?, 1?
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=0; i<mField->getNpx(); i++)
			for(int j=0; j<mField->getNpy()-1; j++)
				HX(i,j) = CHX(i,j)*HX(i,j) 
						- CHXLY(i,j)*( EZX(i,j+1)-EZX(i,j) + EZY(i,j+1)-EZY(i,j));	//Hx�̌v�Z Hx(i, j+1/2) -> Hx[i,j]

				
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=0; i<mField->getNpx()-1; i++)
			for(int j=0; j<mField->getNpy(); j++)
				HY(i,j) = CHY(i,j)*HY(i,j, 0) 
						+ CHYLX(i,j)*( EZX(i+1,j)-EZX(i,j) + EZY(i+1,j)-EZY(i,j) );	//Hy�̌v�Z Hy(i+1/2, j) -> Hy[i,j]

		}
	};

};
#endif //_ST_FDTD_TM_H