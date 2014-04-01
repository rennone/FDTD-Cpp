#ifndef _FDTD_TE_H
#define _FDTD_TE_H
#include "Solver.h"
#include "Object.h"

class FDTD_TE: public Solver{
	typedef Solver super;
protected:
	complex<double> *Ex, *Ey, *Hz, *Hzx, *Hzy;
	double *C_EX, *C_EY, *C_EXLY, *C_EYLX, *C_HZLH;		//LX -> x���� LY -> y����
	double *C_HZX, *C_HZY, *C_HZXLX, *C_HZYLY;
	double *EPS_EX, *EPS_EY, *EPS_HZ;
public:
	FDTD_TE();
	virtual ~FDTD_TE();

	virtual bool calc()=0;
	virtual void draw();
	virtual void field();
	void Initialize();

	void NsScatteredWave(int angle);	//�U���g
	void IncidentWave(int angle);
	void IncidentWaveH(int angle);

	virtual void NTFFindexform(string label, NTFF::output flag = NTFF::REFLEC);

	//�Q�b�^�[
	complex<double>& EX(const int &i, const int &j, const int &t){
		return Ex[pmlIndex(i,j, t)];
	};

	complex<double>& EY(const int &i, const int &j, const int &t){
		return Ey[pmlIndex(i,j, t)];
	};

	complex<double>& HZ(const int &i, const int &j, const int &t){
		return Hz[pmlIndex(i,j, t)];
	};

	complex<double>& EX(const int &i, const int &j){
		return Ex[pmlIndex(i,j)];
	};

	complex<double>& EY(const int &i, const int &j){
		return Ey[pmlIndex(i,j)];
	};

	complex<double>& HZ(const int &i, const int &j){
		return Hz[pmlIndex(i,j)];
	};

	complex<double>& HZX(const int &i, const int &j){
		return Hzx[pmlIndex(i,j)];
	};

	complex<double>& HZY(const int &i, const int &j){
		return Hzy[pmlIndex(i,j)];
	};

	double& CEX(const int &i, const int &j){
		return C_EX[pmlIndex(i,j)];
	};

	double& CEY(const int &i, const int &j){
		return C_EY[pmlIndex(i,j)];
	};

	double& CEXLY(const int &i, const int &j){
		return C_EXLY[pmlIndex(i,j)];
	};

	double& CEYLX(const int &i, const int &j){
		return C_EYLX[pmlIndex(i,j)];
	};

	double& CHZLH(const int &i, const int &j){
		return C_HZLH[pmlIndex(i,j)];
	}

	double& CHZX(const int &i, const int &j){
		return C_HZX[pmlIndex(i,j)];
	}

	double& CHZY(const int &i, const int &j){
		return C_HZY[pmlIndex(i,j)];
	}

	double& CHZXLX(const int &i, const int &j){
		return C_HZXLX[pmlIndex(i,j)];
	}

	double& CHZYLY(const int &i, const int &j){
		return C_HZYLY[pmlIndex(i,j)];
	}

	double& EPSHZ(const int &i, const int &j){
		return EPS_HZ[pmlIndex(i,j)];
	};

	double& EPSEX(const int &i, const int &j){
		return EPS_EX[pmlIndex(i,j)];
	};

	inline double& EPSEY(const int &i, const int &j){
		return EPS_EY[pmlIndex(i,j)];
	};

	void SaveData(string prefix = "");
	void OpenData(string prefix = "");
};
#endif //_FDTD_TE_H