#include "Object.h"
#include "StFDTD_TM.h"

using namespace std;

StFDTD_TM::StFDTD_TM()
:FDTD_TM()
{
	cout << "StFDTD_TM Constructor" << endl;
}

StFDTD_TM::~StFDTD_TM(){
	cout << "StFDTD_TM Destructor" << endl;
}

void StFDTD_TM::field(){
	super::field();
	double sig_x, sig_y, sig_xx, sig_yy;//σx, σx*, σy, σy*
	for(int i=0; i<mField->getNpx(); i++){
		for(int j=0; j<mField->getNpy(); j++){
			sig_x = mField->sigmaX(i,j);
			sig_xx = MU_0_S/EPSILON_0_S * sig_x;
			sig_y = mField->sigmaY(i,j);
			sig_yy = MU_0_S/EPSILON_0_S * sig_y;
			//Δt = 1, μ(i,j) = μ0 
			CEZX(i,j)   = MaxwellCoef(EPSEZ(i,j), sig_x);
			CEZXLX(i,j) = MaxwellCoef2(EPSEZ(i,j), sig_x);

			CEZY(i,j)   = MaxwellCoef(EPSEZ(i,j), sig_y);
			CEZYLY(i,j) = MaxwellCoef2(EPSEZ(i,j), sig_y);

			CHX(i,j)    = MaxwellCoef(MU_0_S, sig_yy);
			CHXLY(i,j)  = MaxwellCoef2(MU_0_S, sig_yy);

			CHY(i,j)    = MaxwellCoef(MU_0_S, sig_xx);
			CHYLX(i,j)  = MaxwellCoef2(MU_0_S, sig_xx);
		}
	}
}

bool StFDTD_TM::calc(){

	CalcE();
	//EZX(mField->getNpx()/2, mField->getNpy()/2) += 0.5*ray_coef*polar(1.0, w_s*time);
	//EZY(mField->getNpx()/2, mField->getNpy()/2) += 0.5*ray_coef*polar(1.0, w_s*time);
	NsScatteredWave(0);
	CalcH();
	if(time > 4000){
		MiePrint(Ez, "Mie_TM2");
		return false;
	}
	
	/*
	if(file){
		file << Ez[index(9*mField->getNx()/10, mField->getNy()/2, +1)] << endl;
	}
	else{		
		file.open("../../Fourie/Fourie/TM_data2.txt");
	}
	*/

	if(time > 7000) return false;
	
	return true;
}

//----------------吸収境界-----------------------//

void StFDTD_TM::absorbing(){
	absorbing_stRL(Ez,0,	LEFT);
	absorbing_stRL(Ez,mField->getNx()-1,	RIGHT);
	absorbing_stTB(Ez,0,	BOTTOM);
	absorbing_stTB(Ez,mField->getNy()-1,	TOP);
}

//----------------周期境界-----------------------//
void StFDTD_TM::cycle(){
	cycle_stRL(Ez,0   ,LEFT);
	cycle_stRL(Ez,mField->getNx()-1,RIGHT);
	//i=0 にi=Nx-2をコピー, i=Nx-1にi=1をコピー
	for(int i=0; i<mField->getNx(); i++){
		EZ(i,0 ,+1) = EZ(i, mField->getNy()-2, +1);
		EZ(i,mField->getNy()-1, +1) = EZ(i,1, +1);
	}
}

/*
	ButtonFactory::setButton("epsx  ", EPSHX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("epsy  ", EPSHY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHYLX  ", CHYLX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHXLY  ", CHXLY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZXLX  ", CEZXLX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZYLY  ", CEZYLY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHY  ", CHY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHX  ", CHX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZX  ", CEZX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZY  ", CEZY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("Hy(right)", HY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2).real());
	ButtonFactory::setButton("Hy(top)", HY(mField->getNpx()/2 ,mField->getNpy()- mField->getNpml()).real());
*/