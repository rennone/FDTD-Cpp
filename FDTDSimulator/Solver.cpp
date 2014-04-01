#include "Solver.h"
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <math.h>
#include <stdio.h> 
#include <string>
#include <algorithm>
#include <glut.h>
#include <fstream>
#include <iostream>

#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define _USE_MATH_DEFINES

Solver::Solver()
	:H_S(1.0), DT_S(1.0)
{
	mField = new Field(2000, 2000, 10, 20); //width, height, ��h, Npml
	LambdaRange    = Range<double>(Nano_S(380), Nano_S(700), Nano_S(5));
	WaveAngleRange = Range<int>   (0, 0, 10);

	SetWaveParameter( LambdaRange.MIN() );
	wave_angle  = WaveAngleRange.MIN();

	time = 0;
	maxStep  = 2500;

	n_s      = new double[mField->getNcel()];	//���ܗ�
	//mModel    = new FazzySlabModel(mField);
	mModel	 = new FazzyMieModel(mField, lambda_s);
	//mModel	 = new FazzyMorphoModel(150, 55, NONSHELF);
	//mModel	 = new FazzyNoModel(mField);
	DataDir		=  "../../../DataSet/";
	WorkingDir  =  "";

	cout << "Solver Constructor" << endl;
}

Solver::~Solver(){
	delete[] n_s;
	delete mModel;
	delete mField;
	cout << "Solver Destructor" << endl;
}

//Bilinear Interpolation���  �����l�̔z��ԍ�������(�l������̔�Œl���o��)
double Solver::bilinear_interpolation(complex<double> *p, double x, double y){
	int i = floor(x);
	int j = floor(y);
	double dx = (x - 1.0*i);	//�����_�ȉ��̒l
	double dy = (y - 1.0*j);	

	return     norm(p[index(i,    j)]) * (1.0-dx)*(1.0-dy)
	         + norm(p[index(i+1,  j)]) * dx*(1.0-dy)
	         + norm(p[index(i,  j+1)]) * (1.0-dx)*dy
	         + norm(p[index(i+1,j+1)]) * dx*dy;
}

//�J���[�}�b�v
Color Solver::color(double phi){
	double range = 2.0;
	Color c;
	double ab_phi = abs(phi);

//	double a = ab_phi < 1 ? (ab_phi <  0.34 ? min(1.0, max(0.0, 3*ab_phi)) : (-1.5*ab_phi+2) ) : 0.5;
	double a = ab_phi < range ? (ab_phi <  range/3.0 ? 3.0/range*ab_phi : (-3.0/4.0/range*ab_phi+1.25) ) : 0.5;
	c.red  = phi > 0 ? a:0;
	c.blue = phi < 0 ? a:0;
	c.green = min(1.0, max(0.0, -3*ab_phi+2));

	return c;
}

void Solver::MiePrint(complex<double>* p, string name){
	printf("Mie print\n");
	ofstream ofs = WriteOpen(name+"Mie");

	if(ofs) {
		for(int i=0; i<=180; i++){
			double _x = 1.2*lambda_s*cos((180-i)*PI/180) + mField->getNx()/2;
			double _y = 1.2*lambda_s*sin((180-i)*PI/180) + mField->getNy()/2;
			double _val = bilinear_interpolation(p,_x,_y);
			ofs << _val << endl;
		}
		printf("output finished\n");
	}
	else{
		printf("No file");
	}
}


//---------------------------------------//
//				 ���˔g				 	 //
//--------------------------------------//
//�_�����̓���
void Solver::pointLightSource(complex<double> *p){
	p[index(mField->getNpx()/2, mField->getNpy()/2, +1)] += ray_coef*polar(1.0, w_s*time);
}

//�������̓���
void Solver::linearLightSource(complex<double> *p){
	for(int i=1; i<mField->getNy()-1; i++)
		p[index(5,i, +1)]  +=ray_coef*polar(1.0, w_s*time);	//todo ���E��ɂ�����Ă����̂�?
}

//�U���g�̌v�Z
void Solver::scatteredWave(complex<double> *p){
	double rad = wave_angle*M_PI/180;	//���W�A���ϊ�
	double _cos = cos(rad), _sin = sin(rad);	//����v�Z����Ǝ��Ԃ����肻��������,������Ă���

	for(int i=mField->getNpml(); i<mField->getNpx(); i++){
		for(int j=mField->getNpml(); j<mField->getNpy(); j++){
			if( N_S(i,j) == 1.0 ) continue;		//���ܗ���1�Ȃ�U���͋N���Ȃ�
			double ikx = k_s*(i*_cos + j*_sin);
			p[index(i,j, +1)] += ray_coef*(1/_pow(N_S(i,j), 2)-1)
				                    *(polar(1.0, ikx-w_s*(time+DT_S))+polar(1.0, ikx-w_s*(time-DT_S))-2.0*polar(1.0, ikx-w_s*time)); 
		}
	}
}


//4�ߖT�𒲂ׂ�
bool Solver::neighber(int _x, int _y){
	for(int i=-2; i<=1; i++)
		if(n_s[index(_min(mField->getNx()-1, _max(0, _x+i%2)), _min(mField->getNy()-1, _max(0,_y+(i+1)%2)) )] != n_s[index(_x,_y)]) return false;

	return true;
}

//-------------------------------------------//
//------------St�z�����E--------------------//
//------------------------------------------//
void Solver::absorbing_stRL(complex<double> *p, int X, enum DIRECT offset){
	double u;
	for(int j=1; j<mField->getNy()-1; j++){
		u = LIGHT_SPEED_S*DT_S/n_s[index(X,j)];
		if(j == 1 || j == mField->getNy()-2)		// �l���̉��͈ꎟ���z�����E
			p[index(X,j, +1)] = p[index(X+offset,j, 0)] + (1- u)/(1+u)*(p[index(X,j, 0)] - p[index(X+offset,j, +1)]);

		else						//����ȊO�͓񎟌��z�����E
			p[index(X,j, +1)] = - p[index(X+offset,j, -1)] 
								     - (1-u)/(1+u)*(p[index(X,j, -1)] + p[index(X+offset,j, +1)]) 
								     +     2/(1+u)*(p[index(X,j,  0)] + p[index(X+offset,j,  0)]) 
								     + u*u/(1+u)/2*( Dy2(p, X,j, 0)   +  Dy2(p, X+offset,j, 0)	);
												        //  dy^2 ��n     +   dy^2 ��b
	}
}

void Solver::absorbing_stTB(complex<double> *p, int Y, enum DIRECT offset){
	double u;
	for(int i=1; i<mField->getNx()-1; i++){
		u = LIGHT_SPEED_S*DT_S/n_s[index(i,Y)];

		if(i==1 || i==mField->getNx()-2)	//�l���̉��͈ꎟ���z�����E
			p[index(i,Y, +1)]    = p[index(i,Y+offset, 0)]    + (1- u)/(1+u)*(p[index(i,Y, 0)]    - p[index(i,Y+offset, +1)]);

		else				//�񎟌��z�����E
			p[index(i,Y, +1)]  = -p[index(i,Y+offset, -1)]    
								  - (1-u)/(1+u)*(p[index(i,Y, -1)] + p[index(i,Y+offset, +1)]) 
							      +     2/(1+u)*(p[index(i,Y, 0)]  + p[index(i,Y+offset, 0)])     
								  + u*u/(1+u)/2*( Dx2(p, i,Y, 0)   + Dx2(p, i,Y+offset, 0) 	);
												  //  dx^2 ��n     +   dx^2 ��b
	}		
}

//-----------------------------------------------------------//
//-------------------------Ns�z�����E------------------------//
//----------------------------------------------------------//
/**���E�̕ǂ�NS�z�����E
** �K�p�z��
** �K�p����x���W
** �E������
*/
void Solver::absorbing_nsRL(complex<double> *p, int X, enum DIRECT offset){
	double kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	double ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;
	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int j=1; j<mField->getNy()-1; j++){
		u1 = tan(w_b/n_s[index(X,j)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/N_S(X,j)), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(j == 1 || j == mField->getNy()-2)		// �l���̉��͈ꎟ���z�����E
			p[index(X,j, +1)] = p[index(X+offset,j, 0)] + (1- u1)/(1+u1)*(p[index(X,j, 0)] - p[index(X+offset,j, +1)]);

		else						//����ȊO�͓񎟌��z�����E
			p[index(X,j, +1)] = - p[index(X+offset,j, -1)] 
								- (1-u1)/(1+u1)*(p[index(X,j, -1)] + p[index(X+offset,j, +1)]) 
								+     2/(1+u1)*(p[index(X,j,  0)] + p[index(X+offset,j,  0)]) 
								+ u2*u2/(1+u1)/2*( Dy2(p, X,j, 0)   +  Dy2(p, X+offset,j, 0)	);
												        //  dy^2 ��n     +   dy^2 ��b
	}
}

/**�㉺�̕ǂ�NS�z�����E
** �K�p�z��
** �K�p����y���W
** �ォ����
*/
void Solver::absorbing_nsTB(complex<double> *p, int Y, enum DIRECT offset){
	double kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	double ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;
	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int i=1; i<mField->getNx()-1; i++){
		u1 = tan(w_b/n_s[index(i,Y)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/n_s[index(i,Y)]), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(i==1 || i==mField->getNx()-2)	//�l���̉��͈ꎟ���z�����E
			p[index(i,Y, +1)]    = p[index(i,Y+offset, 0)]    + (1- u1)/(1+u1)*(p[index(i,Y, 0)] - p[index(i,Y+offset, +1)]);

		else				//�񎟌��z�����E
			p[index(i,Y, +1)]  = -p[index(i,Y+offset, -1)]    
								  - (1-u1)/(1+u1)*(p[index(i,Y, -1)] + p[index(i,Y+offset, +1)]) 
							      +     2/(1+u1)*(p[index(i,Y, 0)]  + p[index(i,Y+offset, 0)])     
								  + u2*u2/(1+u1)/2*( Dx2(p, i,Y, 0)   + Dx2(p, i,Y+offset, 0) 	);
												  //  dx^2 ��n     +   dx^2 ��b
	}		
}

//----------------------------------------------//
//----------------�������E-----------------------//
//-----------------------------------------------//

void Solver::cycle_stRL(complex<double> *p, int X, enum DIRECT offset){
	double u;
	for(int i=1; i<mField->getNy()-1; i++){
		u = LIGHT_SPEED_S*DT_S/n_s[index(X,i)];
			//�񎟌��z�����E
		p[index(X,i, +1)] = - p[index(X+offset,i, -1)]     
							- (1-u)/(1+u)*(p[index(X,i, -1)] + p[index(X+offset,i, +1)])    
							+     2/(1+u)*(p[index(X,i,  0)] + p[index(X+offset,i,  0)]) 
							+ u*u/(1+u)/2*(  Dy2(p,X,i, 0)   +   Dy2(p,X+offset,i, 0) );
												//  dy^2 ��n     +   dy^2 ��b
	}		
}

//---------------------------�`��------------------------------//
//������(0,0), �E�オ(Nx,Ny)
void Solver::draw(Complex *p, Complex *q){
	double N = max(mField->getNx(),mField->getNy());
	double ws = 2.0*MAIN_WINDOW_X/WINDOW_W/N;
	double hs = 2.0*MAIN_WINDOW_H/WINDOW_H/N;
	for (int i = mField->getNpml(); i < mField->getNpx()-mField->getNpml(); i++){
		for (int j = mField->getNpml(); j < mField->getNpy()-mField->getNpml(); j++){
			int x = i-mField->getNpml();
			int y = j-mField->getNpml();
			Color c = color( norm(p[index(i,j)] + q[index(i,j)]) );
			//Color c = color(30.0*(p[index(i,j)].real() + q[index(i,j)].real()));
			glColor3d(c.red, c.green, c.blue);
			if(j==mField->getNpy()/2 || i==mField->getNpx()/2) glColor3d(1,1,1);
			
			glRectd(x*ws-1, y*hs-1, (x+1.0)*ws-1, (y+1.0)*hs-1);
				
		}
	}

	//draw_model();
}

//---------------------------�`��------------------------------//
//������(0,0), �E�オ(Nx,Ny)
void Solver::draw(complex<double> *p){
	double N = max(mField->getNx(),mField->getNy());
	double ws = 2.0*MAIN_WINDOW_X/WINDOW_W/N;
	double hs = 2.0*MAIN_WINDOW_H/WINDOW_H/N;
	for (int i = mField->getNpml(); i < mField->getNpx()-mField->getNpml(); i++){
		for (int j = mField->getNpml(); j < mField->getNpy()-mField->getNpml(); j++){
			int x = i-mField->getNpml();
			int y = j-mField->getNpml();
			//Color c = color( norm(p[index(i,j)]) );
			Color c = color(30.0*p[index(i,j)].real());
			glColor3d(c.red, c.green, c.blue);
			if(j==mField->getNpy()/2 || i==mField->getNpx()/2) glColor3d(1,1,1);
			
			glRectd(x*ws-1, y*hs-1, (x+1.0)*ws-1, (y+1.0)*hs-1);
				
		}
	}

	draw_model();
}

//�U���̂̕`��
void Solver::draw_model(){
	double N = max(mField->getNx(),mField->getNy());
	double ws = 2.0*MAIN_WINDOW_X/WINDOW_W/N;
	double hs = 2.0*MAIN_WINDOW_H/WINDOW_H/N;
	for (int i = mField->getNpml(); i < mField->getNpx()-mField->getNpml(); i++){
		for (int j = mField->getNpml(); j < mField->getNpy()-mField->getNpml(); j++){
			int x = i-mField->getNpml();
			int y = j-mField->getNpml();
			//�}�����E
			const double n = N_S(i,j);	//������,���ܗ������������Ă͂����Ȃ�

			if(n == 1.0) continue;	//���ܗ���1�Ȃ�Ƃ΂�
			glColor4d(0.7/n, 0.7/n, 0.7/n, 0.1);
			glRectd(x*ws-1, y*hs-1, (x+1.0)*ws-1, (y+1.0)*hs-1);	
		}
	}
}

//-----------------�f�[�^�̕ۑ�----------------------//
void Solver::save_data(complex<double> *data, string name){
	ofstream out = WriteOpen(name);

	for(int k=-1; k<2; k++)
		for(int i=0; i<mField->getNx(); i++)
			for(int j=0; j<mField->getNy(); j++)
				out << data[index(i,j, k)] << endl;

}

void Solver::open_data(complex<double> *data, string name){
	ifstream in = ReadOpen(name);

	for(int k=-1; k<2; k++)
		for(int i=0; i<mField->getNx(); i++)
			for(int j=0; j<mField->getNy(); j++){
				if(!in) {throw  "cannot open file" ;}
				in >> data[index(i,j, k)];
			}
}



//---------------------------------------------//
//--------------PML----------------------------//
//--------------------------------------------//

/*
Color Solver::color2(double a){
	Color c;
	c.red  = max(0.0, min(1.0, a-1.0));
	c.blue = max(0.0, min(-a+1, 1.0));
	c.green = max(0.0, min(1.0, -2*abs(a-1)+1));
	return c;
}

*/