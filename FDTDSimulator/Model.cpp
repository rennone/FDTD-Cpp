#include"Model.h"
#include"Field.h"

/*---------------------------------------------*/
/*--------------�~��Mie�U��--------------------*/
/*---------------------------------------------*/
FazzyMieModel::FazzyMieModel(Field *f, double _r):
FazzyModel(f),r(_r)
{
	ep = 1.6*1.6*EPSILON_0_S;
		cout << r << endl;
}

string FazzyMieModel::mkdir(string root){
	_mkdir((root + "Mie").c_str());

	string name = "Mie/" + to_s((int)(mField->cellToNano(r))) +"nm,"+ mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f){

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)�����_�ɃV�t�g
	double _y = my - 0.5*mField->getNy();

	//���S�Ƃ̋��������a+��2/2�Z���ȏ�Ȃ�, ���S�ɔ}���̊O(�O�̂���, ���a+1 �ȏォ���ׂĂ���)
	if(_x*_x + _y*_y >= pow(r+1, 2.0))

		return EPSILON_0_S;

	//���S�Ƃ̋�����, ���a-��2/2�Z���ȓ��Ȃ�, ���S�ɔ}���̊O
	if(_x*_x + _y*_y <= pow(r-1, 2.0))
		return ep;

	double s=0;

	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1)
		for(double j=-16+0.5; j<16; j+=1)
			if(pow(_x+a*i/32.0, 2.0) + pow(_y+b*j/32.0, 2.0) <= r*r)
				s+=1; 
	s /= 32.0*32.0;
	return ep*s + EPSILON_0_S*(1-s);
}

/*---------------------------------------------*/
/*--------------���w��-------------------------*/
/*---------------------------------------------*/
FazzySlabModel::FazzySlabModel(Field* f):
FazzyModel(f), ep1(2.0*2.0*EPSILON_0_S), ep2(EPSILON_0_S), width1(250), width2(50)
{
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f){
//��100nm����,250nm�Ԋu��50nm�̃X���u�����Ă���
//���w��
	
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();

	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	int k    = (int)(mField->cellToNano(mx) - 100)%250;
	double l =      (mField->cellToNano(mx) - 100)/250;

	if( k > 0 && k <=50 && l < 5)
		return ep1;
	else
		return ep2;

}

string FazzySlabModel::mkdir(string root){
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

/*---------------------------------------------*/
/*--------------�����t�H��--------------------*/
/*---------------------------------------------*/
FazzyMorphoModel::FazzyMorphoModel(Field* f, double _h0, double _h1, enum STRUCTURE kind):
FazzyModel(f), shelf(kind)
{
	num = 8;				//�ςݏd�˂鐔
	// ep[1] = 3.5*3.5*EPSILON_0_S;	//�U�d��
	// ep[0] = 1.45*1.45*EPSILON_0_S;
	ep[1] = 1.56*1.56*EPSILON_0_S;	//�U�d��
	ep[0] = 1.0*1.0*EPSILON_0_S;
	width = mField->nanoToCell(150);	//������150�ŌŒ�

	min = mField->nanoToCell(120);
	max = mField->nanoToCell(120);
	height[1] = min;
	height[0] = min;
	cout << min << endl;
	cout << max << endl;
}

double FazzyMorphoModel::calcEPS(const double &x, const double &y, enum INTEG f){
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	//N_X/2�𒆐S����,�����`�����Ɍ݂��Ⴂ�ɔz�u
	int dis = height[0] + height[1];
	int oy = (mField->getNy() - num*dis)/2.0;
	int ox = mField->getNx()/2.0;
	double _x = mx-ox;	//ox,oy�����W�̌��_��
	double _y = my-oy;
	
	//���f���̍��E���㉺�ɗ���Ă���Δ}���̊O
	if(abs(_x) > width+1 || abs(_y-1.0*num*dis/2.0) > num*dis/2.0+1)
		return EPSILON_0_S;

	double s[2]={0,0};
	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1){
		for(double j=-16+0.5; j<16; j+=1){
			double sx = _x + a*i/32.0;
			double sy = _y + b*j/32.0;
			if( abs(sx) > width  || abs(sy-1.0*num*dis/2.0) > 1.0*num*dis/2.0) continue;

			bool k = ((int)sy%dis+0.5) > height[0];	//���E��Ŕ�ׂȂ��悤��0.5�����Ă���(0���傫��1�����Ȃ牽�ł�����)
			//bool k =  (floor(sy/ dis)*dis < sy) && ( sy < floor(sy/ dis)*dis + height[0]);

			if (sx < 0 && shelf)
				k = !k;		//���E�Ŕ��], �݂��Ⴂ�łȂ������甽�]���Ȃ�

			s[k] +=1;
		}
	}
	s[0] /= 32.0*32.0;
	s[1] /= 32.0*32.0;
	return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

string FazzyMorphoModel::mkdir(string root){
	string label = "Morpho(" + to_s(sqrt(ep[0])) + "," + to_s(sqrt(ep[1])) + ")M=" + to_s(num);
	//  string label = "Morpho";
	_mkdir((root + label).c_str());
	string name;

	if(shelf)
		name = label + "/" + to_s((int)(mField->cellToNano(height[0]))) + "nm" + mField->getStringCellInfo() ;
	else
		name = label + "/" + to_s((int)(mField->cellToNano(height[0]))) + "nm(nonShelf)" + mField->getStringCellInfo();

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}


/*--------------------------------*/
/*-----------���f���Ȃ�-----------*/
/*--------------------------------*/
bool FazzyMorphoModel::update(int dh){
	height[0] += (int) mField->nanoToCell(dh);
	height[1] += (int) mField->nanoToCell(dh);

	if(height[0] > max)
		return false;

	return true;
}

string FazzyNoModel::mkdir(string root){
	_mkdir((root + "NoModel").c_str());

	string name = "NoModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}