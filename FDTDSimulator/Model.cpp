#include"Model.h"
#include"Field.h"

/*---------------------------------------------*/
/*--------------円柱Mie散乱--------------------*/
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
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f){

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	if(_x*_x + _y*_y >= pow(r+1, 2.0))

		return EPSILON_0_S;

	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
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
/*--------------多層膜-------------------------*/
/*---------------------------------------------*/
FazzySlabModel::FazzySlabModel(Field* f):
FazzyModel(f), ep1(2.0*2.0*EPSILON_0_S), ep2(EPSILON_0_S), width1(250), width2(50)
{
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f){
//左100nmから,250nm間隔で50nmのスラブを入れていく
//多層膜
	
	double mx = x - mField->getNpml(); //計算領域内へ写像
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
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

/*---------------------------------------------*/
/*--------------モルフォ蝶--------------------*/
/*---------------------------------------------*/
FazzyMorphoModel::FazzyMorphoModel(Field* f, double _h0, double _h1, enum STRUCTURE kind):
FazzyModel(f), shelf(kind)
{
	num = 8;				//積み重ねる数
	// ep[1] = 3.5*3.5*EPSILON_0_S;	//誘電率
	// ep[0] = 1.45*1.45*EPSILON_0_S;
	ep[1] = 1.56*1.56*EPSILON_0_S;	//誘電率
	ep[0] = 1.0*1.0*EPSILON_0_S;
	width = mField->nanoToCell(150);	//横幅は150で固定

	min = mField->nanoToCell(120);
	max = mField->nanoToCell(120);
	height[1] = min;
	height[0] = min;
	cout << min << endl;
	cout << max << endl;
}

double FazzyMorphoModel::calcEPS(const double &x, const double &y, enum INTEG f){
	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	//N_X/2を中心軸に,長方形を横に互い違いに配置
	int dis = height[0] + height[1];
	int oy = (mField->getNy() - num*dis)/2.0;
	int ox = mField->getNx()/2.0;
	double _x = mx-ox;	//ox,oyを座標の原点に
	double _y = my-oy;
	
	//モデルの左右か上下に離れていれば媒質の外
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

			bool k = ((int)sy%dis+0.5) > height[0];	//境界上で比べないように0.5足している(0より大きく1未満なら何でもいい)
			//bool k =  (floor(sy/ dis)*dis < sy) && ( sy < floor(sy/ dis)*dis + height[0]);

			if (sx < 0 && shelf)
				k = !k;		//左右で反転, 互い違いでなかったら反転しない

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

	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}


/*--------------------------------*/
/*-----------モデルなし-----------*/
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
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}