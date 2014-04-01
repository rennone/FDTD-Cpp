#define _USE_MATH_DEFINES
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include "Simulator.h"
#include "EventState.h"
#include "Solver.h"

using namespace std;
Simulator::Simulator(){
	//Field *mField = new Field(2000, 2000, 10, 20); //width, height, Δh, Npml
	//Model *mModel	 = new FazzyMieModel(mField, lambda_s);
	solv = new StFDTD_TM();
	solv->field();
}

Simulator::~Simulator(){ 
	cout << "Simulator Destructor" << endl;
	ButtonFactory::deleteAllButton();	//ボタンの削除, ウィンドウを閉じて終了した際に,すでに解放されたボタンを解放する恐れがある.(この順番だと大丈夫かも)
	delete solv;						//solverの削除
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
	solv->draw();					//シミュレーション状況描画
	ButtonFactory::draw();		//ボタンを描画
}
