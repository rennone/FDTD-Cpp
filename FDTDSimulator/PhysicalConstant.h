#ifndef _PHYSICAL_H
#define _PHYSICAL_H

#include <math.h>
/*-----------------------------*/
/*�����萔�̃V�~�����[�V�����l */
/*-----------------------------*/
const double LIGHT_SPEED_S  = 0.7;	//����
const double EPSILON_0_S    = 1.0;	//�^��̗U�d��
const double MU_0_S	        = 1.0/0.7/0.7;	//�^��̓������@//Debug���[�h����,0.7��Light_SPEED_S�Ə�����,0�ŏ���������Ă��܂�.(Release�@�Ȃ�ok)
const double Z_0_S	        = sqrt(1.0/0.7/0.7/1.0);	//�g���C���s�[�_���X z = ��(��/��)
#endif