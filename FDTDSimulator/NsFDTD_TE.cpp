#include"NsFDTD_TE.h"
#include <omp.h>

NsFDTD_TE::NsFDTD_TE()
:FDTD_TE()
{
	cout << "NsFDTD_TE Constructor" << endl;
};

NsFDTD_TE::~NsFDTD_TE(){
	cout << "NsFDTD_TE Destructor" << endl;
};

bool NsFDTD_TE::calc(){	

	CalcE();	//電界の計算

	//NsScatteredWave(wave_angle);
	absorbing();

	CalcH();		//磁界の計算 Hz(i+1/2, j+1/2) -> Hz[i,j]

	pointLightSource(Hz);

	if(time > maxStep)
		return EndTask();
	
	return true;
};


bool NsFDTD_TE::EndTask(){
	cout << "End Task" << endl;
	string label = "";

	NTFFindexform(label, NTFF::NTFFDATA | NTFF::TOTAL);

	//終了条件の確認
	if( !Terminate())
		return false;

	ReStart();

	return true;
}

void NsFDTD_TE::field(){		
	super::field();	//誘電率の設定
	setWorkingDirPass(mModel->mkdir(getDataDirPass()));	//データを保存するディレクトリの設定

	R_M = NsCoef();		//差分演算用の計算定数の設定
	R_P = 1.0-R_M;

	//σ=0を用いて最適化した定数の計算
	double mu;
	double sig = 0;	//真空の誘電率, 導電率, 透磁率
	for(int i=0; i<mField->getNx(); i++){
		for(int j=0; j<mField->getNy(); j++){
			mu = MU_0_S;
			double u_hz = sin(w_s / sqrt(EPSHZ(i,j) / EPSILON_0_S) * DT_S / 2) / sin(k_s*DT_S/2);
			double u_ex = sin(w_s / sqrt(EPSEX(i,j) / EPSILON_0_S) * DT_S / 2) / sin(k_s*DT_S/2);
			double u_ey = sin(w_s / sqrt(EPSEY(i,j) / EPSILON_0_S) * DT_S / 2) / sin(k_s*DT_S/2);
			CEX(i,j)   = 1.0;
			CEXLY(i,j) = u_ex*sqrt(mu/EPSEX(i,j));

			CEY(i,j)   = 1.0;
			CEYLX(i,j) = u_ey*sqrt(mu/EPSEY(i,j));	

			CHZLH(i,j) = u_hz*sqrt(EPSHZ(i,j)/mu);;

		}
	}
	cout << "field" << endl;
}

void NsFDTD_TE::absorbing(){
	//H(i,j)だったやつ（今までちゃんと動いていたやつ）
	
	absorbing_nsRL(Ey, 0,	 LEFT);
	absorbing_nsRL(Ey, mField->getNx()-2, RIGHT);
	absorbing_nsTB(Ex, 0,    BOTTOM);
	absorbing_nsTB(Ex, mField->getNy()-2, TOP);
	/*
	absorbing_nsRL(Ey, 1,	 LEFT);
	absorbing_nsRL(Ey, mField->getNx()-1, RIGHT);
	absorbing_nsTB(Ex, 1,    BOTTOM);
	absorbing_nsTB(Ex, mField->getNy()-1, TOP);
	*/
}



//double a = 0; //α = σ/(2ε) より, 今はσ=0としているのでαも0
//double u = sqrt( (_pow(sin(sqrt(w_s*w_s - a*a)*DT_S/2 ),2)+_pow( sinh(a*DT_S/2),2) )/ (_pow(sin(k_s*DT_S/2),2)*cosh(a*DT_S))  );
// a = 0, tanh(0) = sinh(0) = 0, cosh(0) = 1　を用いて最適化する


