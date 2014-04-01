#define _USE_MATH_DEFINES
#include "Object.h"
#include "FDTD_TE.h"
#include "Vec.h"
#include <math.h> 
using namespace std;

FDTD_TE::FDTD_TE()
:Solver()
{
	cout << "FDTD_TE Constructor" << endl;
	Hz  = new complex<double>[mField->getNcel()];		//Hz(i+0.5, j+0.5) �� Hz(i,j)���Ӗ����� �̈�m��
	Ex  = new complex<double>[mField->getNcel()];		//Ex(i+0.5,j) �� Ex(i,j)���Ӗ�����
	Ey  = new complex<double>[mField->getNcel()];		//Ey(i,j+0.5) �� Ey(i,j)���Ӗ�����
	Hzx= new complex<double>[mField->getNcel()];		//Ey(i,j+0.5) �� Ey(i,j)���Ӗ�����
	Hzy= new complex<double>[mField->getNcel()];		//Ey(i,j+0.5) �� Ey(i,j)���Ӗ�����

//�v�Z�p�萔�z��
	C_HZLH  = new double[mField->getNcel()];	//CHzlh(i+0.5, j+0.5) �� CHZLH(i,j)
	C_EX    = new double[mField->getNcel()];	//CEx(i+0.5, j) �� CEX(i,j)
	C_EY    = new double[mField->getNcel()];	//CEy(i, j+0.5) �� CEY(i,j)
	C_EXLY  = new double[mField->getNcel()];
	C_EYLX  = new double[mField->getNcel()];
	C_HZX   = new double[mField->getNcel()];
	C_HZY   = new double[mField->getNcel()];
	C_HZXLX = new double[mField->getNcel()];
	C_HZYLY = new double[mField->getNcel()];
	EPS_HZ  = new double[mField->getNcel()];
	EPS_EX  = new double[mField->getNcel()];
	EPS_EY  = new double[mField->getNcel()];


	//������
	for(int i=0; i<mField->getNcel(); i++)
		Ex[i] = Ey[i] = Hz[i] = Hzx[i] = Hzy[i] = 0;

	//todo �萔�͏��������Ȃ��Ă���?
	// C_HZX[i] = C_HZY[i] = C_HZXLX[i] = C_HZYLY[i] = C_HZLH[i] = C_EX[i] = C_EY[i] = C_EXLY[i] = C_EYLX[i] = 0;
}

FDTD_TE::~FDTD_TE(){
	delete[] Hz;
	delete[] Ex;
	delete[] Ey;
	delete[] Hzx;
	delete[] Hzy;

	delete[] C_HZLH;
	delete[] C_EX;
	delete[] C_EY;
	delete[] C_EXLY;
	delete[] C_EYLX;
	delete[] C_HZX;
	delete[] C_HZY;
	delete[] C_HZXLX;
	delete[] C_HZYLY;

	delete[] EPS_HZ;
	delete[] EPS_EX;
	delete[] EPS_EY;

	cout << "FDTD_TE Destructor" << endl;
}


void FDTD_TE::draw(){
	super::draw(Ey);
}

void FDTD_TE::OpenData(string prefix){
	string str = prefix+ "data/" + getWaveData();	//�f�[�^�I�[�v��
	cout << str << endl;

	open_data(Hz, str +"Hz");
	open_data(Ex, str +"Ex");
	open_data(Ey, str +"Ey");
}

void FDTD_TE::SaveData(string prefix){
	string str = MakeDir(prefix+ "data");
	cout << str << endl;
	save_data(Hz, str +  "/" + getWaveData() + "Hz");
	save_data(Ex, str + "/" + getWaveData()  + "Ex");
	save_data(Ey, str + "/" + getWaveData()  + "Ey");
}

void FDTD_TE::Initialize(){
	super::Initialize();
	//�̈揉����
	for(int i=0; i<3*mField->getNcel(); i++)
			Ex[i] = Ey[i] = Hz[i] = 0;
}

void FDTD_TE::field(){
	setWorkingDirPass(MakeDir("TE"));

	for(int i=0; i<mField->getNpx(); i++){
		for(int j=0; j<mField->getNpy(); j++){
			EPSEX(i,j) = mModel->calcEPS(i+0.5,j    , D_Y);
			EPSEY(i,j) = mModel->calcEPS(i    ,j+0.5, D_X);
			EPSHZ(i,j) = 0.5*(mModel->calcEPS(i+0.5,j+0.5, D_X) + mModel->calcEPS(i+0.5,j+0.5, D_Y));
			N_S(i,j)   = sqrt(mModel->calcEPS(i+0.5,j+0.5) / EPSILON_0_S );		//���ܗ� n = ��(��/ ��0) (�����ł� �� = ��0�Ƃ��Ă���̂�)
		}
	}
}

//�U���g
void FDTD_TE::NsScatteredWave(int ang){
	double rad = ang*M_PI/180;	//���W�A���ϊ�
	double _cos = cos(rad), _sin = sin(rad);	//����v�Z����Ǝ��Ԃ����肻��������,������Ă���
	//double a = (1-exp(-_pow(0.01*time,2)));		//�s�A���ɓ��˂����̂�h���ׂ̌W��
	double a = ray_coef;
	double n, u0, u1, _n;
	for(int i=0; i<mField->getNx(); i++){
		for(int j=0; j<mField->getNy(); j++){
			if(EPSEY(i,j) == EPSILON_0_S && EPSEX(i,j) == EPSILON_0_S) continue;
			double ikx = k_s*(i*_cos + j*_sin);

			n  = sqrt(EPSEY(i,j)/EPSILON_0_S);				//todo �� = 0���O��̘b
			u0 = sin(w_s/n*H_S/2.0)  / sin(k_s*H_S/2.0);
			u1 = sin(w_s/n*DT_S/2.0) / sin(k_s*n*H_S/2.0);
			_n = u0/u1;
			EY(i,j, +1) += a*_cos*(1.0/(_n*n) - 1.0)*(polar(1.0, ikx - w_s*(time+DT_S)) - polar(1.0, ikx - w_s*time));

			n  = sqrt(EPSEX(i,j)/EPSILON_0_S);				//todo �� = 0���O��̘b
			u0 = sin(w_s/n*H_S/2.0)  / sin(k_s*H_S/2.0);
			u1 = sin(w_s/n*DT_S/2.0) / sin(k_s*n*H_S/2.0);
			_n = u0/u1;
			EX(i,j, +1) += a*_sin*(1.0/(_n*n) - 1.0)*(polar(1.0, ikx - w_s*(time+DT_S)) - polar(1.0, ikx - w_s*time));

		}
	}
}

//���˔g
void FDTD_TE::IncidentWave(int ang){
	double rad = ang*M_PI/180;	//���W�A���ϊ�
	double _cos = cos(rad), _sin = sin(rad);	//����v�Z����Ǝ��Ԃ����肻��������,������Ă���
	for(int i=0; i<mField->getNx(); i++){
		for(int j=0; j<mField->getNy(); j++){
			double ikx = k_s*(i*_cos + j*_sin);
			EY(i,j, +1) += _cos*polar(1.0, ikx - w_s*(time+DT_S));
			EX(i,j, +1) += _sin*polar(1.0, ikx - w_s*(time+DT_S));
		}
	}
}

void FDTD_TE::IncidentWaveH(int ang){
	double rad = ang*M_PI/180;	//���W�A���ϊ�
	double _cos = cos(rad), _sin = sin(rad);	//����v�Z����Ǝ��Ԃ����肻��������,������Ă���
	for(int i=0; i<mField->getNx(); i++){
		for(int j=0; j<mField->getNy(); j++){
			double ikx = k_s*(i*_cos + j*_sin);
			HZ(i,j, +1) += k_s/(w_s*MU_0_S) * polar(1.0, ikx - w_s*(time+DT_S))*(_sin - _cos);
		}
	}
}


void FDTD_TE::NTFFindexform(string label, NTFF::output flag){
	cout << "NTF start" << endl;

	int cx = mField->getNx()/2;
	int cy = mField->getNy()/2;

	double r0 = 1.0e6;
    complex<double> iu( 0.0, 1.0 );											// imaginary unit
    complex<double> Coef = sqrt( iu*k_s/(8*M_PI*r0) ) * exp( iu*k_s*r0 );	// common coefficient

    int offset = 5;		    // closed line offset
	int lt,rt, tp, bm;		//�Ȗʂ̏ꏊ
	tp = mField->getNy()-offset;			//�ォ��-5
	bm = offset;			//������5
	rt = mField->getNx() - offset;		//�E����-5
	lt  = offset;			//������5

	double sum = 0;
	const int max_angle = 360;	//�ǂ̊p�x�܂ŕ��z�����߂邩, 180��360
	  double strength[max_angle];
    for ( int phi=0; phi<max_angle; phi++ ) {
        double rad = phi * M_PI/180.0;
		rad = M_PI - rad;							//todo �Ȃ���TE��180�����炷
        Vec2<double> r( cos(rad), sin(rad) );	        // eye vector
          
        complex<double> Nx( 0, 0 );
        complex<double> Ny( 0, 0 );
        complex<double> Lz( 0, 0 );
        Vec2<double> r2;
        complex<double> C_EX, C_EY, C_HZ;

        // (left,bottom) -> (right,bottom)
        for ( int i=lt; i<rt; i++ ) {
            r2    = Vec2<double>( i-cx + 0.5, bm-cy);	//���S����Z���܂ł̋���
            C_HZ  = 0.5*( HZ(i,bm, +1) + HZ(i,bm-1, +1) );
            C_EX  = EX(i,bm, +1);
            Nx   -= C_HZ * exp( iu * k_s * In_prod(r,r2) );
            Lz   -= C_EX * exp( iu * k_s * In_prod(r,r2) );
        }

        // (right,bottom) -> (right,top)
        for ( int j=bm; j<tp; j++ ) {
			r2   = Vec2<double>( rt-cx, j-cy + 0.5 );     
            C_HZ  = 0.5*( HZ(rt,j, +1) + HZ(rt+1,j, +1) );  
            C_EY  = 0.25 * ( EY(rt,j, +1) + EY(rt-1,j, +1) + EY(rt,j+1, +1) + EY(rt-1,j+1, +1));  
            Ny -= C_HZ * exp( iu * k_s * In_prod(r,r2) );
            Lz -= C_EY * exp( iu * k_s * In_prod(r,r2) );
        }

        // (right,top) -> (left,top)
        for ( int i=lt; i<rt; i++ ) {
            r2    = Vec2<double>( i-cx+0.5, tp-cy );
            C_HZ  = 0.5*(HZ(i,tp, +1) + HZ(i,tp+1, +1));
            C_EX  = 0.25 * ( EX(i,tp, +1) + EX(i,tp-1, +1) + EX(i+1,tp, +1) + EX(i+1,tp-1, +1));
            Nx += C_HZ * exp( iu * k_s * In_prod(r,r2) );
            Lz += C_EX * exp( iu * k_s * In_prod(r,r2) );
        }            

        // (left,top) -> (left,bottom)
        for ( int j=bm; j<tp; j++ ) {
            r2    = Vec2<double>( lt-cx, j-cy+0.5 );         
            C_HZ  = 0.5*(HZ(lt,j, +1) + HZ(lt-1,j, +1));
            C_EY  = 0.25 * ( EY(lt,j, +1) + EY(lt-1,j, +1) + EY(lt,j+1, +1) + EY(lt-1,j+1, +1));  
            Ny += C_HZ * exp( iu * k_s * In_prod(r,r2) );
            Lz += C_EY * exp( iu * k_s * In_prod(r,r2) );
        }
            
        // Get Ephi
        complex<double> Ephi = Coef * ( -Z_0_S*( -Nx*sin(rad) + Ny*cos(rad) ) + Lz );
          
		strength[phi] = norm(Ephi);
		sum += norm(Ephi);
    }

	//NTFF�̌��ʂ��o��
	if( (flag & NTFF::NTFFDATA) == NTFF::NTFFDATA ){
		ofstream fp = WriteOpen(MakeDir("NTFF") + label + getWaveData());
		for(int i=0; i<max_angle; i++)
			fp << strength[i] << endl;
	}

	//NTFF�̌��ʂ̑��a���o��
	if( (flag & NTFF::TOTAL) == NTFF::TOTAL){
		ofstream fp= WriteOpen(MakeDir("NTFF") + label + "WaveAngleStrength", DATAFILE::ADD); //false�͒ǋL���[�h
		fp << "(" <<(int)Inv_Nano_S(lambda_s) << "," << wave_angle << ")  " << sum << endl;		//�g�����Ƃ̑��a���o��
	}

	//���˗����o��
	if( (flag & NTFF::REFLEC) == NTFF::REFLEC){
		ofstream fp = WriteOpen(MakeDir("Reflection") + label + getWaveData());
		for(int i = 0; i < max_angle; i++)
			fp << strength[i] / sum << endl;
	}
	cout << "NTF finish" << endl;
}
