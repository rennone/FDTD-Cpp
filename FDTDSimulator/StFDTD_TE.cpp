#include"StFDTD_TE.h"
#include "Object.h"

using namespace std;

StFDTD_TE::StFDTD_TE()
:FDTD_TE()
{
	cout << "StFDTD_TE Constructor" << endl;
}

StFDTD_TE::~StFDTD_TE(){
	cout << "StFD_TE Destructor" << endl;
}

bool StFDTD_TE::calc(){
	CalcE();
	pointLightSource(Ey);
	CalcH();
	if(time > 3000){
		MiePrint(Ey, "Mie_TE2");
		return false;
	}

	return true;
}

void StFDTD_TE::field(){
	super::field();

	double sig_x, sig_xx, sig_y, sig_yy; 	//��x, ��x*, ��y, ��y*
	for(int i=0; i<mField->getNpx();i++){
		for(int j=0; j<mField->getNpy();j++){
			sig_x = mField->sigmaX(i,j);
			sig_xx = MU_0_S/EPSILON_0_S * sig_x;
			sig_y = mField->sigmaY(i,j);
			sig_yy = MU_0_S/EPSILON_0_S * sig_y;

			//��t = 1, ��(i,j) = ��0 
			CEX(i,j)  =  MaxwellCoef(EPSEX(i,j), sig_y);
			CEXLY(i,j)= MaxwellCoef2(EPSEX(i,j), sig_y);

			CEY(i,j)  = MaxwellCoef(EPSEY(i,j), sig_x);
			CEYLX(i,j)= MaxwellCoef2(EPSEY(i,j), sig_x);

			CHZX(i,j)  = MaxwellCoef(MU_0_S, sig_xx);
			CHZXLX(i,j)= MaxwellCoef2(MU_0_S, sig_xx);

			CHZY(i,j)  = MaxwellCoef(MU_0_S, sig_yy);
			CHZYLY(i,j)= MaxwellCoef2(MU_0_S, sig_yy);
		}
	}
}

//----------------Mur�̋z�����E-----------------------//
void StFDTD_TE::absorbing(){
	absorbing_stTB(Ex, mField->getNy()-2, TOP);
	absorbing_stTB(Ex, 0,    BOTTOM);
	absorbing_stRL(Ey, mField->getNx()-2, RIGHT);
	absorbing_stRL(Ey, 0,    LEFT);
}

//----------------�������E��������-----------------------//
void StFDTD_TE::cycle(){
	
	cycle_stRL(Hz,0   ,LEFT);
	cycle_stRL(Hz,mField->getNx()-1,RIGHT);
	//i=0 ��i=mField->getNx()-2���R�s�[, i=mField->getNx()-1��i=1���R�s�[
	for(int i=0; i<mField->getNx(); i++){
		HZ(i,0 ,+1)	   = HZ(i, mField->getNy()-2, +1);
		HZ(i,mField->getNy()-1, +1) = HZ(i,1, +1);
	}
}

/*
	ButtonFactory::setButton("sigx  ", EPSEX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("sigy  ", EPSEY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEYLX  ", CEYLX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEXLY  ", CEXLY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHZXLX  ", CHZXLX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHZYLY  ", CHZYLY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEY  ", CEY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEX  ", CEX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHZX  ", CHZX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHZY  ", CHZY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("Ey(right)", EY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2).real());
	ButtonFactory::setButton("Ey(top)", EY(mField->getNpx()/2 ,mField->getNpy()- mField->getNpml()).real());
*/