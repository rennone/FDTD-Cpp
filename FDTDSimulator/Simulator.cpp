#define _USE_MATH_DEFINES
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include "Simulator.h"
#include "EventState.h"
#include "Solver.h"

using namespace std;
Simulator::Simulator(){
	//Field *mField = new Field(2000, 2000, 10, 20); //width, height, ��h, Npml
	//Model *mModel	 = new FazzyMieModel(mField, lambda_s);
	solv = new StFDTD_TM();
	solv->field();
}

Simulator::~Simulator(){ 
	cout << "Simulator Destructor" << endl;
	ButtonFactory::deleteAllButton();	//�{�^���̍폜, �E�B���h�E����ďI�������ۂ�,���łɉ�����ꂽ�{�^����������鋰�ꂪ����.(���̏��Ԃ��Ƒ��v����)
	delete solv;						//solver�̍폜
};

int Simulator::calc()
{
	try{
		solv->nextTime();
		return solv->calc();
	}
	catch(char *err){
		cout << err << endl;
		return 0;
	}
	return 1;
}

void Simulator::draw()
{
//	return;
	if( ((int)solv->getTime()) %20 != 0) return;
	solv->draw();					//�V�~�����[�V�����󋵕`��
	ButtonFactory::draw();		//�{�^����`��
}
