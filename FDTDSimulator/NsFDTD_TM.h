#ifndef _NS_FDTD_TM_
#define _NS_FDTD_TM_
#include "FDTD_TM.h"

class NsFDTD_TM: public FDTD_TM{
	typedef FDTD_TM super;
	double R_P, R_M;
public:
	NsFDTD_TM();
	~NsFDTD_TM();
	bool calc();
	void field();
private:

	complex<double> Dx2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i,j+1, t)] + p[index(i, j-1, t)] - p[index(i-1,j+1, t)] - p[index(i-1, j-1, t)])/2.0;
	};

	complex<double> Dy2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i+1,j, t)] + p[index(i-1, j, t)] - p[index(i+1,j-1, t)] - p[index(i-1, j-1, t)])/2.0;
	}

	//todo ���E�ߖT��S-FDTD���g���ĂȂ�
	void CalcE(){
	#ifdef _OPENMP
	#pragma omp parallel for
	#endif
		//�d�E�̌v�Z
		for(int i=1; i<mField->getNx()-1; i++){
			for(int j=1; j<mField->getNy()-1; j++){

				EZ(i,j,+1) = CEZ(i,j)*EZ(i,j, 0)
					+ CEZLH(i,j)*(
					R_P*( HY(i,j, 0) - HY(i-1,j, 0) - ( HX(i,j, 0) - HX(i,j-1, 0) )  )	//dx(Hy) - dy(Hy)
					+R_M*(  Dx2_n(Hy,i,j, 0) - Dy2_n(Hx,i,j, 0) )
					);
			}
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
		for(int i=0; i<mField->getNx(); i++){
			for(int j=0; j<mField->getNy()-1; j++){
				HX(i,j, +1) = HX(i,j, 0) - CHXLY(i,j)*( EZ(i,j+1, +1) - EZ(i,j, +1) );	//Hx�̌v�Z Hx(i, j+1/2) -> Hx[i,j]
			}
		}
				
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=0; i<mField->getNx()-1; i++){
			for(int j=0; j<mField->getNy(); j++){
				HY(i,j, +1) = HY(i,j, 0) + CHYLX(i,j)*( EZ(i+1,j, +1) - EZ(i,j, +1) );	//Hy�̌v�Z Hy(i+1/2, j) -> Hy[i,j]
			}
		}
			
		}
	};

	void absorbing();	//�z�����E
	void PML_E();
	void PML_H();

	bool EndTask();		//1��̃V�~�����[�V�������I������Ƃ��̏���

	void ReStart(){
		super::Initialize();
		field();
	}
	void printDebugCalcE(int,int,int,int);
	void printDebugCalcHx(int,int,int,int);
	void printDebugCalcHy(int,int,int,int);
};

//�v�����g�f�o�b�O
inline void NsFDTD_TM::printDebugCalcE(int i, int j, int _i, int _j){
	if(i != _i || j!=_j) return;
	cout << "===========pringDebug CalcE TM Mode===============" << endl; 
	cout << "(" + to_s(i) + "," + to_s(j) + ")" << endl;
	cout << "EZ(i,j,+1) = " << EZ(i,j,+1) << endl;
	cout << "EZ(i,j, 0) = " << EZ(i,j,0) << endl;
	cout << "EZ(i,j,-1) = " << EZ(i,j,-1) << endl;
	cout << "CE(i,j) = " << CEZ(i,j) << endl;
	cout << "CEZLH(i,j) = " <<  CEZLH(i,j) << endl;
	cout << "HY(i  ,j  , 0) = " << HY(i,j,0) << endl;
	cout << "HY(i-1,j  , 0) = " << HY(i-1,j,0) << endl;
	cout << "HY(i  ,j+1, 0) = " << HY(i,j+1,0) << endl;
	cout << "HY(i  ,j-1, 0) = " << HY(i,j-1,0) << endl;
	cout << "HY(i-1,j+1, 0) = " << HY(i-1,j+1,0) << endl;
	cout << "HY(i-1,j-1, 0) = " << HY(i-1,j-1,0) << endl;

	cout << "HX(i,j,0) = " << HX(i,j,0) << endl;

	cout << "HX(i  ,j-1,0) = " << HY(i,j-1,0) << endl;
	cout << "HX(i+1,j  ,0) = " << HX(i+1,j,0) << endl;
	cout << "HX(i-1,j  ,0) = " << HY(i-1,j,0) << endl;
	cout << "HX(i+1,j-1,0) = " << HX(i+1,j-1,0) << endl;
	cout << "HX(i-1,j-1,0) = " << HY(i-1,j-1,0) << endl;
	cout << "R_P = " << R_P << endl;
	cout << "R_M = " << R_M << endl;
	cout << "==================================================" << endl; 
};

//�v�����g�f�o�b�O
inline void NsFDTD_TM::printDebugCalcHx(int i, int j, int _i, int _j){
	if(i != _i || j!=_j) return;
	cout << "===========pringDebug CalcHx TM Mode===============" << endl; 
	cout << "(" + to_s(i) + "," + to_s(j) + ")" << endl;
	cout << "HXi,j,+1) = "  << HX(i,j,+1) << endl;
	cout << "HX(i,j, 0) = " << HX(i,j,0) << endl;
	cout << "HX(i,j,-1) = " << HX(i,j,-1) << endl;

	cout << "CHXLY(i,j) = " << CHXLY(i,j) << endl;
	cout << "EZ(i  ,j  , 0) = " << EZ(i,j,0) << endl;
	cout << "EZ(i  ,j+1, 0) = " << EZ(i,j+1,0) << endl;
	cout << "==================================================" << endl; 
}

//�v�����g�f�o�b�O
inline void NsFDTD_TM::printDebugCalcHy(int i, int j, int _i, int _j){
	if(i != _i || j!=_j) return;
	cout << "===========pringDebug CalcHy TM Mode===============" << endl; 
	cout << "(" + to_s(i) + "," + to_s(j) + ")" << endl;
	cout << "HYi,j,+1) = "  << HY(i,j,+1) << endl;
	cout << "HY(i,j, 0) = " << HY(i,j,0) << endl;
	cout << "HY(i,j,-1) = " << HY(i,j,-1) << endl;

	cout << "CHYLX(i,j) = " << CHYLX(i,j) << endl;
	cout << "EZ(i  ,j  , 0) = " << EZ(i,j,0) << endl;
	cout << "EZ(i+1,j  , 0) = " << EZ(i+1,j,0) << endl;
	cout << "==================================================" << endl; 
}


#endif

/*
				EZ(i,j,+1) = CEZ(i,j)*EZ(i,j, 0)
					+ CEZLH(i,j)*(
					R_P*( HY(i+1,j, 0) - HY(i,j, 0) - ( HX(i,j+1, 0) - HX(i,j, 0) )  )	//dx(Hy) - dy(Hy)
					+R_M*(  Dx2_n(Hy,i,j, 0) - Dy2_n(Hx,i,j, 0) )
					);
					*/
//HX(i,j, +1) = HX(i,j, 0) - CHXLY(i,j)*( EZ(i,j, +1) - EZ(i,j-1, +1) );	//Hx�̌v�Z Hx(i, j+1/2) -> Hx[i,j]
			
//HY(i,j, +1) = HY(i,j, 0) + CHYLX(i,j)*( EZ(i,j, +1) - EZ(i-1,j, +1) );	//Hy�̌v�Z Hy(i+1/2, j) -> Hy[i,j]