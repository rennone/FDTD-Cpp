#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include"NsFDTD_TM.h"
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

//---------------�R���X�g���N�^, �f�X�g���N�^ -----------------//
NsFDTD_TM::NsFDTD_TM()
:FDTD_TM()
{
	cout << "NsFDTD_TM Constructor" << endl;
};

NsFDTD_TM::~NsFDTD_TM(){
	cout << "NsFDTD_TM Destructor" << endl;
};


bool NsFDTD_TM::calc(){
	CalcE();	//�d��̌v�Z
	//NsScatteredWave(wave_angle);	//�U���g�̓���
	pointLightSource(Ez);
	absorbing();					//�z�����E

	CalcH();	//����̌v�Z
	ButtonFactory::setButton("EZ", norm(EZ(mField->getNx()/2+10, mField->getNy()/2)));
	if(time>maxStep)
		return EndTask();

	return true;
};

//----------------�I�����̎d��------------------------//
bool NsFDTD_TM::EndTask(){
	cout << "End Task" << endl;

//	string label = "";
//	NTFFindexform("", NTFF::NTFFDATA | NTFF::TOTAL );	// label -> "" �ɂ������Ǔ������m�F���ĂȂ�.
	
	//�I�������̊m�F
	if( !Terminate() )
		return false;

	ReStart();
	return true;
}

//--------------�v�Z�W���̐ݒ�---------------------//
void NsFDTD_TM::field(){		
	super::field();										//�U�d���̐ݒ�
	//setWorkingDirPass(mModel->mkdir(getDataDirPass()));	//�f�[�^��ۑ�����f�B���N�g���̐ݒ�
	R_M = NsCoef();		//�������Z�p�̌v�Z�萔�̐ݒ�
	R_P = 1.0-R_M;
	for(int i=0; i<mField->getNx(); i++){
		for(int j=0; j<mField->getNy(); j++){
			double mu = MU_0_S;
			
			double u_ez = sin(w_s/sqrt(EPSEZ(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
			double u_hx = sin(w_s/sqrt(EPSHX(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
			double u_hy = sin(w_s/sqrt(EPSHY(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
			
			
			CEZ(i,j)   = 1;							//(1 - tanh(a*h_s)) / (1 + tanh(a*H_S));
			CEZLH(i,j) = u_ez*sqrt(mu/EPSEZ(i,j));		//u' ��(��/��)*(1/(1+tanh(a*H_S)))				
			CHXLY(i,j) = u_hx*sqrt(EPSHX(i,j)/mu);		//u' ��(��/��)
			CHYLX(i,j) = u_hy*sqrt(EPSHY(i,j)/mu);		//u' ��(��/��)
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

/* Field�̌W���ݒ�̕���
	double u_ez = sin(w_s/sqrt(EPSEZ(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
	double u_hx = sin(w_s/sqrt(EPSHX(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
	double u_hy = sin(w_s/sqrt(EPSHY(i,j)/EPSILON_0_S)*DT_S/2)/ sin(k_s*DT_S/2);
*/
/*
//MorphoModel���������,�����I�ɔz�u����Ă��邩�̊m�F�p, field��SetEPS�̌�ɔz�u
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
	//�v�Z�p�萔�̐ݒ� �ԈႦ�Ă����ۂ�
	double kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	double ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;

*/

	/*
	double sig = 0;	//���d��=0�ōl����
	double a = 0;	//�� = ��/(2��) ���, ���̓�=0�Ƃ��Ă���̂Ń���0
	//double u = sqrt( (_pow(sin(sqrt(w_s*w_s - a*a)*DT_S/2 ),2)+_pow( sinh(a*DT_S/2),2) )/ (_pow(sin(k_s*DT_S/2),2)*cosh(a*DT_S))  );

	// a = 0, tanh(0) = sinh(0) = 0, cosh(0) = 1�@��p���čœK������
	//double u = sin(w_s*DT_S/2)/ sin(k_s*DT_S/2);
	*/

/*
	if(time <= maxStep){
		OpenData(name);
		time = maxStep+1;
	}
	else{
		//SaveData(name);		//�V�~�����[�g�����f�[�^�͕ۑ�
	}
*/