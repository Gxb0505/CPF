#include "CPF.hpp"

// ����
int		Constant::_HSVResolution	= 10;
int		Constant::_Resolution		= 110;
double	Constant::_HSVRatio			= 0.8;
int		Constant::_particleNum		= 200;

// ��ʼ����������
double	Constant::_priorNoiseWidth	= 4;
double	Constant::_priorNoiseHeight	= 4;
double	Constant::_priorNoiseScale	= 0.1;

Rand	gRand;