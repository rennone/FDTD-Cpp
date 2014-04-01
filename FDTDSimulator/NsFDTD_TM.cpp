#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include"NsFDTD_TM.h"
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

//---------------コンストラクタ, デストラクタ -----------------//
NsFDTD_TM::NsFDTD_TM()
:FDTD_TM()
{
	cout << "NsFDTD_TM Constructor" << endl;
};

NsFDTD_TM::~NsFDTD_TM(){
	cout << "NsFDTD_TM Destructor" << endl;
};


bool NsFDTD_TM::calc(){
	CalcE();	//電場の計算
	//NsScatteredWave(wave_angle);	//散乱波の入射
	pointLightSource(Ez);
	absorbing();					//吸収境界

	CalcH();	//磁場の計算
	ButtonFactory::setButton("EZ", norm(EZ(mField->getNx()/2+10, mField->getNy()/2)));
	if(time>maxStep)
		return EndTask();

	return true;
};

//----------------終了時の仕事------------------------//
bool NsFDTD_TM::EndTask(){
	cout << "End Task" << endl;

//	string label = "";
//	NTFFindexform("", NTFF::NTFFDATA | NTFF::TOTAL );	// label -> "" にしたけど動くか確認してない.
	
	//終了条件の確認
	if( !Terminate() )
		return false;

	ReStart();
	return true;
}

//--------------計算係数の設定---------------------//
void NsFDTD_TM::field(){		
	super::field();										//誘電率の設定
	//setWorkingDirPass(mModel->mkdir(getDataDirPass()));	//データを保存するディレクトリの設定
	R_M = NsCoef();		//差分演算用の計算定数の設定
	R_P = 1.0-R_M;
	for(int i=0; i<mField->getNx(); i++){
		for(int j=0; j<mField->getNy(); j++){
			double mu = MU_0_S;
			
			double u_ez = sin(w_s/sqrt(EPSEZ(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
			double u_hx = sin(w_s/sqrt(EPSHX(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
			double u_hy = sin(w_s/sqrt(EPSHY(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
			
			
			CEZ(i,j)   = 1;							//(1 - tanh(a*h_s)) / (1 + tanh(a*H_S));
			CEZLH(i,j) = u_ez*sqrt(mu/EPSEZ(i,j));		//u' √(μ/ε)*(1/(1+tanh(a*H_S)))				
			CHXLY(i,j) = u_hx*sqrt(EPSHX(i,j)/mu);		//u' √(ε/μ)
			CHYLX(i,j) = u_hy*sqrt(EPSHY(i,j)/mu);		//u' √(ε/μ)
		}
	}
	cout << "Field" << endl;
}

void NsFDTD_TM::absorbing(){
	absorbing_nsRL(Ez, 0,	 LEFT);
	absorbing_nsRL(Ez, mField->getNx()-1, RIGHT);
	absorbing_nsTB(Ez, 0,    BOTTOM);
	absorbing_nsTB(Ez, mField->getNy()-1, TOP);
}

/* Fieldの係数設定の部分
	double u_ez = sin(w_s/sqrt(EPSEZ(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
	double u_hx = sin(w_s/sqrt(EPSHX(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
	double u_hy = sin(w_s/sqrt(EPSHY(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
*/
/*
//MorphoModelがきちんと,周期的に配置されているかの確認用, fieldのSetEPSの後に配置
	cout << "confirm" << endl;
	int s=0;
	double bf = 1.0;
	for(int i=0; i < mField->getNy(); i++){
		if(N_S(mField->getNx()/2, i) == bf){
			s++;
		}
		else{
			cout << bf << "=" << s << endl;
			s = 1;
		}
		bf = N_S(mField->getNx()/2, i);
	}
	cout << "confirm end" << endl;
	*/

/*
	//計算用定数の設定 間違えてたっぽい
	double kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	double ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;

*/

	/*
	double sig = 0;	//導電率=0で考える
	double a = 0;	//α = σ/(2ε) より, 今はσ=0としているのでαも0
	//double u = sqrt( (_pow(sin(sqrt(w_s*w_s - a*a)*DT_S/2 ),2)+_pow( sinh(a*DT_S/2),2) )/ (_pow(sin(k_s*DT_S/2),2)*cosh(a*DT_S))  );

	// a = 0, tanh(0) = sinh(0) = 0, cosh(0) = 1　を用いて最適化する
	//double u = sin(w_s*DT_S/2)/ sin(k_s*DT_S/2);
	*/

/*
	if(time <= maxStep){
		OpenData(name);
		time = maxStep+1;
	}
	else{
		//SaveData(name);		//シミュレートしたデータは保存
	}
*/