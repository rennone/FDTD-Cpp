#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <glut.h>
#include "Simulator.h"
#include "EventState.h"
#include<Windows.h>
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

//----------------------------------------------------
// �ϐ��̐錾
//----------------------------------------------------
const int WindowPositionX = 200;  //��������E�B���h�E�ʒu��X���W
const int WindowPositionY = 200;  //��������E�B���h�E�ʒu��Y���W
const char WindowTitle[] = "�V�~�����[�V����";  //�E�B���h�E�̃^�C�g��

Simulator *wave_simulator = new Simulator();

//----------------------------------------------------
// �֐��v���g�^�C�v�i��ɌĂяo���֐����ƈ����̐錾�j
//----------------------------------------------------
void Initialize(void);   //�����ݒ莞�ɌĂяo���֐�
void Idle(void);         //�A�C�h�����ɌĂяo���֐�
void Display(void);      //��ʕ`�掞�ɌĂяo���֐�
void Mouse(int _button , int _state , int _x , int _y);
void Motion(int _x,int _y);
void KeyBoard(unsigned char key, int x, int y);
void SpecialKeyBoard(int key, int x, int y);

//----------------------------------------------------
// ���C���֐�
//----------------------------------------------------

int main(int argc, char *argv[]){
//	cout.setf(ios::fixed, ios::floatfield);	//�Œ菬���_, 18�����x�Ŏw��
//	cout.precision(18);
	
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );	//���������[�N���o�p
	glutInit(&argc, argv);                                     //���̏�����
	glutInitWindowPosition(WindowPositionX, WindowPositionY);  //�E�B���h�E�̈ʒu�̎w��
	glutInitWindowSize(WINDOW_W, WINDOW_H);					//�E�B���h�E�T�C�Y�̎w��
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE); //�f�B�X�v���C���[�h�̎w��
	glutCreateWindow(WindowTitle);                             //�E�B���h�E�̍쐬
	glutIdleFunc(Idle);                                        //�v���O�����A�C�h����Ԏ��ɌĂяo�����֐�
	glutDisplayFunc(Display);                                  //�`�掞�ɌĂяo�����֐����w�肷��i�֐����FDisplay�j
	glutMouseFunc(Mouse);
	glutKeyboardFunc(KeyBoard);
	glutSpecialFunc(SpecialKeyBoard);
	glutMotionFunc(Motion);
	Initialize();                                              //�����ݒ�̊֐����Ăяo��
	glutMainLoop();  
	return 0;

}
//----------------------------------------------------
// �����ݒ�̊֐�
//----------------------------------------------------
void Initialize(void){

}
//----------------------------------------------------
// �A�C�h�����ɌĂяo�����֐�
//----------------------------------------------------
void Idle(){

	//if( getSpecialKeyState() == GLUT_KEY_UP){
		if(!wave_simulator->calc()){
			delete wave_simulator;
			exit(0);
		}
	//}

	if( getKeyState() == _KEY_ESC){
		delete wave_simulator;
		exit(0);
	}
	
	Sleep(1);
  setKeyState(-1);			//�L�[�{�[�h������
  setSpecialKeyState(-1);
  glutPostRedisplay(); //glutDisplayFunc()���P����s����
}

void Mouse(int _button , int _state , int _x , int _y){
	setMouseState( _button,  _state, _x, _y);
}

void Motion(int _x,int _y){
	setPos(_x,_y);
}

void KeyBoard(unsigned char key, int x, int y){
	setKeyState(key);
}

void SpecialKeyBoard(int key, int x, int y){
	setSpecialKeyState(key);
}
//----------------------------------------------------
// �`��̊֐�
//----------------------------------------------------
void Display(void) {
	glClear(GL_COLOR_BUFFER_BIT);

	wave_simulator->draw();

	glutSwapBuffers(); //glutInitDisplayMode(GLUT_DOUBLE)�Ń_�u���o�b�t�@�����O�𗘗p��
	glFlush();
}
