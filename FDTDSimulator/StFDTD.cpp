#include "FDTD.h"

StFDTD::StFDTD()
:FDTD()
{
	cout << "StFDTD Constructor" << endl;

	ButtonFactory::addButton("time", 0, UNCHECK);
	ButtonFactory::addButton("phi" , 0, UNCHECK);
	ButtonFactory::addButton("str" , 0, UNCHECK);

}

StFDTD::~StFDTD(){
	cout << "StFDTD Destructor" << endl;
}

//�v�Z�p�萔�̐ݒ�
void StFDTD::field(){
	//Mie_Slub();
	for(int i=0; i<mField->getNx(); i++)
		for(int j=0; j<mField->getNy(); j++)
			np[index(i,j)] = _pow(LIGHT_SPEED_S*DT_S/H_S/n_s[index(i,j)], 2);	
}

//phi���v�Z
bool StFDTD::calc(){
	time += DT_S;		//���̃X�e�b�v�̌v�Z

	//���E�ȊO
	for(int i=1; i<mField->getNx()-1; i++)
		for(int j=1; j<mField->getNy()-1; j++)
			phi[index(i,j, +1)] = np[index(i,j)]*( Dx2(phi, i,j, 0) + Dy2(phi, i,j, 0) ) - phi[index(i,j, -1)] + 2.0*phi[index(i,j,0)];

	//absorbing();	//�z�����E
	cycle();	//�������E

	linearLightSource(phi);			//i=5�̈ʒu��, ���Ăɔz�u
	//phi[index(mField->getNx()/2, mField->getNy()/2, +1)] += 5*(1-exp(-0.01*time*time))*sin(- w_s*time);

	//�{�^���̒l�X�V
	ButtonFactory::setButton("time", time);
	ButtonFactory::setButton("phi", phi[index(9*mField->getNx()/10, mField->getNy()/2, +1)].real());
	ButtonFactory::setButton("str", norm(phi[index(9*mField->getNx()/10, mField->getNy()/2)]));

	if((int)time%1000 == 0) cout << "time=" << time << endl;	

	/*
	if(file)
		file << phi[index(9*mField->getNx()/10, mField->getNy()/2, +1)] << endl;	//�L�^����
	else
		file.open("../../Fourie/Fourie/data.txt");
		*/
	/*
	if(time > 2000){
		MiePrint(phi, "Mie_phi.txt");
		return false;
	}
	*/


	if(time > 7000) return false;
	

	return true;
}



//�z�����E
void StFDTD::absorbing(){
	absorbing_stRL(phi,0,	 LEFT);
	absorbing_stRL(phi,mField->getNx()-1, RIGHT);
	absorbing_stTB(phi,0,    BOTTOM);
	absorbing_stTB(phi,mField->getNy()-1, TOP);	
}


//�������E
void StFDTD::cycle(){
	cycle_stRL(phi, 0, LEFT);
	cycle_stRL(phi, mField->getNx()-1, RIGHT);
	//i=0 ��i=mField->getNx()-2���R�s�[, i=mField->getNx()-1��i=1���R�s�[
	for(int i=0; i<mField->getNx(); i++){
		phi[index(i,0 ,+1)] = phi[index(i, mField->getNy()-2, +1)];
		phi[index(i,mField->getNy()-1, +1)] = phi[index(i,1, +1)];
	}

}
