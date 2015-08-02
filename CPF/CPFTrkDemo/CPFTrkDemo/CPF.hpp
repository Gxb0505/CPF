

#ifndef __CPF__HPP__INCLUDE__
#define __CPF__HPP__INCLUDE__

/************************************************************************/
/* STL库                                                                */
/************************************************************************/
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
using	namespace	std;

/************************************************************************/
/* OpenCV库                                                             */
/************************************************************************/
#include <opencv2\opencv.hpp> 
using namespace cv; 

/************************************************************************/
/* 其他C库                                                              */
/************************************************************************/
#include <stdarg.h>
#include <time.h>
#include <math.h>

/************************************************************************/
/* 常用函数的定义                                                       */
/************************************************************************/
// 常数控制
class	Constant
{
public:
	static	int				_HSVResolution;					// HSV直方图一个长度 10
	static	int				_Resolution;					// HSV直方图长度 110=10*10+10
	static	double			_HSVRatio;						// HSV直方图HS和V的比率
	static	int				_particleNum;					// 目标粒子个数
	static	double			_priorNoiseWidth;				// 先验噪声
	static	double			_priorNoiseHeight;				// 先验噪声
	static	double			_priorNoiseScale;				// 先验噪声
};

// 全局函数
class	GlobalMethod{
public:
#define			__MATH__PI__			3.1415926
	/************************************************************************
	GblRectCheckBound		检测矩形边界
	GblSmRate				计算相似率
	GblCrossRate			计算交叉
	GblContain				包含 含两个矩形的最小矩形
	GblCross				相交 两个矩形的相交的矩形
	GblGetCenter			得到矩形中心
	************************************************************************/
	static void				GblRectCheckBound		( CvRect* ioRect , int imgWidth , int imgHeight ){
		if (ioRect->x <0)							ioRect->x		= 0;
		if (ioRect->x >imgWidth-1)					ioRect->x		= imgWidth-1;
		if (ioRect->x+ioRect->width <0)				ioRect->width	= 0;
		if (ioRect->x+ioRect->width >imgWidth-1)	ioRect->width	= imgWidth-1 - ioRect->x;

		if (ioRect->y <0)							ioRect->y		= 0;
		if (ioRect->y >imgHeight-1)					ioRect->y		= imgHeight-1;
		if (ioRect->y+ioRect->height <0)			ioRect->height	= 0;
		if (ioRect->y+ioRect->height >imgHeight-1)	ioRect->height	= imgHeight-1 - ioRect->y;
	}

	static	float			GblSmRate				( CvRect& rectFir , CvRect& rectSec ){
		float smWidth , smHeight;
		smWidth	= ( float )min( rectFir.width , rectSec.width ) / ( max( rectFir.width , rectSec.width ) + 0.000001 );
		smHeight= ( float )min( rectFir.height, rectSec.height) / ( max( rectFir.height, rectSec.height) + 0.000001 );
		return sqrt( smWidth * smHeight );
	};

	static	float			GblCrossRate			( CvRect& rectFir , CvRect& rectSec ){
		CvRect	cross	= GblCross( rectFir , rectSec );
		CvRect	contain	= GblContain( rectFir , rectSec );
		return	( float ) (cross.width*cross.height) / (contain.width*contain.height);
	}

	static	CvRect			GblContain				( CvRect& rectFir , CvRect& rectSec ){
		CvRect rectRes;
		rectRes.x		= min( rectFir.x , rectSec.x );
		rectRes.y		= min( rectFir.y , rectSec.y );
		rectRes.width	= max( rectFir.x + rectFir.width , rectSec.x + rectSec.width ) - rectRes.x;
		rectRes.height	= max( rectFir.y + rectFir.height, rectSec.y + rectSec.height) - rectRes.y;
		return rectRes;
	}

	static	CvRect			GblCross				( CvRect& rectFir , CvRect& rectSec ){
		CvRect rectRes;
		rectRes.x = 0;	rectRes.y = 0;	rectRes.width = 0;	rectRes.height = 0;
		if( rectFir.x > ( rectSec.x + rectSec.width ) ){
			return rectRes;
		} else if( rectSec.x > ( rectFir.x + rectFir.width ) ){
			return rectRes;
		} else if( rectFir.y > ( rectSec.y + rectSec.height ) ){
			return rectRes;
		} else if( rectSec.y > ( rectFir.y + rectFir.height ) ){
			return rectRes;
		} else {
			rectRes.x		= max( rectFir.x , rectSec.x );
			rectRes.y		= max( rectFir.y , rectSec.y );
			rectRes.width	= min( rectFir.x + rectFir.width , rectSec.x + rectSec.width ) - rectRes.x;
			rectRes.height	= min( rectFir.y + rectFir.height, rectSec.y + rectSec.height) - rectRes.y;
			return rectRes;
		}
	}

	static	CvPoint			GblGetCenter			( CvRect& rect ){
		return cvPoint( rect.x + rect.width/2 , rect.y + rect.height/2 );
	}
	
	static	CvRect			GblLocalRect			( const IplImage* pImg , CvRect rect , double x , double y , double scale ){
		CvRect rstRect;
		rstRect.x		= cvFloor( x - rect.width / 2 );
		rstRect.y		= cvFloor( y - rect.height / 2 );
		rstRect.width	= cvFloor( rect.width*scale );
		rstRect.height	= cvFloor( rect.height*scale );
		GlobalMethod::GblRectCheckBound( &rstRect , pImg->width , pImg->height );
		return rstRect;
	}

};


// 产生随机数的函数
class	Rand{
#ifndef		__MOD__MAX__
#define		__MOD__MAX__		10000
#endif
public:
	Rand()
	{
		srand ( time(NULL) );
	};
public:
	double				GaussRand			( void ){
		static double V2, fac;
		static int phase = 0;
		double S, Z, U1, U2, V1;
		if( phase ){
			Z = V2 * fac;
		} else {
			do {
				U1 = (double)rand() / RAND_MAX;
				U2 = (double)rand() / RAND_MAX;
				V1 = 2 * U1 - 1;
				V2 = 2 * U2 - 1;
				S = V1 * V1 + V2 * V2;
			} while(S >= 1);
			fac = sqrt (-2 * log(S) / S);
			Z = V1 * fac;
		}
		phase = 1 - phase;
		return Z;
	}
	double				GeneralRand			( void ){
		return (double)( rand() % __MOD__MAX__ ) / __MOD__MAX__;
	}
	double				GeneralRandP5		( void ){
		return (double)( rand() % __MOD__MAX__ ) / __MOD__MAX__ - 0.5;
	}
	int					GeneralRand			( int mod ){
		return rand() % mod;
	}
};
extern	Rand	gRand;

// 计算直方图函数
class	CalHist
{
public:
	/************************************************************************/
	/* 计算HSV直方图                                                        */
	/************************************************************************/
	static void			CalHSVHist			( const IplImage* pImg 
											, const CvRect& rect 
											, double* hist )
	{
		double*	hsHist = new( double[Constant::_HSVResolution*Constant::_HSVResolution] );
		double*	vHist = new( double[Constant::_HSVResolution] );
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hsHist[dh] = 0;
		}
		for( int dh=0; dh<Constant::_HSVResolution; ++dh ){
			vHist[dh] = 0;
		}
		double		pixCount = .0;
		CvScalar	scaPixel;
		int			binIdxH, binIdxS, binIdxV;
		for( int xi=rect.x; xi<rect.x+rect.width; xi++ ){
			for( int yi=rect.y; yi<rect.y+rect.height; yi++ ){
				scaPixel = cvGet2D( pImg , yi , xi );			
				binIdxH	 = cvFloor(scaPixel.val[0]*Constant::_HSVResolution/256);
				binIdxS	 = cvFloor(scaPixel.val[1]*Constant::_HSVResolution/256);
				binIdxV	 = cvFloor(scaPixel.val[2]*Constant::_HSVResolution/256);
				hsHist[binIdxH*Constant::_HSVResolution+binIdxS] += 1;
				vHist[binIdxV] += 1;
				pixCount++;
			}
		}
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hsHist[dh] = hsHist[dh]/pixCount;
		}
		for( int dh=0; dh<Constant::_HSVResolution; ++dh ){
			vHist[dh] = vHist[dh]/pixCount;
		}
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hist[dh] = hsHist[dh] * Constant::_HSVRatio;
		}
		int	vStartInd = Constant::_HSVResolution*Constant::_HSVResolution;
		for( int dv=0; dv<Constant::_HSVResolution; ++dv ){
			hist[dv+vStartInd] = vHist[dv] * (1-Constant::_HSVRatio);
		}
		delete	hsHist;
		delete	vHist;
	}

	/************************************************************************/
	/* Hist Model (包括背景)                                                */
	/************************************************************************/
	static void			CalHSVHist			( const IplImage* pImg 
											, const CvRect& rect 
											, double* hist 
											, double* bgHist )
	{
		double*	hsHist = new( double[Constant::_HSVResolution*Constant::_HSVResolution] );
		double*	vHist = new( double[Constant::_HSVResolution] );
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hsHist[dh] = 0;
		}
		for( int dh=0; dh<Constant::_HSVResolution; ++dh ){
			vHist[dh] = 0;
		}
		double		pixCount = .0;
		CvScalar	scaPixel;
		int			binIdxH, binIdxS, binIdxV;
		for( int xi=rect.x; xi<rect.x+rect.width; xi++ ){
			for( int yi=rect.y; yi<rect.y+rect.height; yi++ ){
				scaPixel = cvGet2D( pImg , yi , xi );			
				binIdxH	 = cvFloor(scaPixel.val[0]*Constant::_HSVResolution/256);
				binIdxS	 = cvFloor(scaPixel.val[1]*Constant::_HSVResolution/256);
				binIdxV	 = cvFloor(scaPixel.val[2]*Constant::_HSVResolution/256);
				hsHist[binIdxH*Constant::_HSVResolution+binIdxS] += 1;
				vHist[binIdxV] += 1;
				pixCount++;
			}
		}
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hsHist[dh] = hsHist[dh]/pixCount;
		}
		for( int dh=0; dh<Constant::_HSVResolution; ++dh ){
			vHist[dh] = vHist[dh]/pixCount;
		}
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hist[dh] = hsHist[dh] * Constant::_HSVRatio;
		}
		int	vStartInd = Constant::_HSVResolution*Constant::_HSVResolution;
		for( int dv=0; dv<Constant::_HSVResolution; ++dv ){
			hist[dv+vStartInd] = vHist[dv] * (1-Constant::_HSVRatio);
		}
		
		// 估计背景直方图
		CvRect	externRect;
		externRect.x	  = rect.x - rect.width / 2;
		externRect.y	  = rect.y - rect.height / 2;
		externRect.width  = rect.width * 2;
		externRect.height = rect.height* 2;
		GlobalMethod::GblRectCheckBound( &externRect , pImg->width , pImg->height );

		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hsHist[dh] = 0;
		}
		for( int dh=0; dh<Constant::_HSVResolution; ++dh ){
			vHist[dh] = 0;
		}
		pixCount = .0;
		for( int xi=externRect.x; xi<rect.x+externRect.width; xi++ ){
			for( int yi=externRect.y; yi<rect.y+externRect.height; yi++ ){
				if( xi>=rect.x && yi>=rect.y && xi<=(rect.x+rect.width) && yi<=(rect.y+rect.height) ){
					
				} else {
					scaPixel = cvGet2D( pImg , yi , xi );			
					binIdxH	 = cvFloor(scaPixel.val[0]*Constant::_HSVResolution/256);
					binIdxS	 = cvFloor(scaPixel.val[1]*Constant::_HSVResolution/256);
					binIdxV	 = cvFloor(scaPixel.val[2]*Constant::_HSVResolution/256);
					hsHist[binIdxH*Constant::_HSVResolution+binIdxS] += 1;
					vHist[binIdxV] += 1;
					pixCount++;
				}
			}
		}
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hsHist[dh] = hsHist[dh]/pixCount;
		}
		for( int dh=0; dh<Constant::_HSVResolution; ++dh ){
			vHist[dh] = vHist[dh]/pixCount;
		}
		for( int dh=0; dh<Constant::_HSVResolution*Constant::_HSVResolution; ++dh ){
			hist[dh] = hsHist[dh] * Constant::_HSVRatio;
		}
		vStartInd = Constant::_HSVResolution*Constant::_HSVResolution;
		for( int dv=0; dv<Constant::_HSVResolution; ++dv ){
			bgHist[dv+vStartInd] = vHist[dv] * (1-Constant::_HSVRatio);
		}

		delete	hsHist;
		delete	vHist;
	}

	static double		CalHistCoffDist		( const double* h1, const double* h2 ){
		double	dis = 0.0;
		for( int dh=0; dh<Constant::_Resolution; ++dh ){
			dis += sqrt( h1[dh] * h2[dh] );
		}
		dis = 1 - dis;
		return	dis;
	}

	static double		CalHistCoffDist	( const double* hc, const double* ho , const double* hb ){
		double	disCO = 0.0 , disCB = 0.0 , dis = 0.0;
		for( int dh=0; dh<Constant::_Resolution; ++dh ){
			disCO += sqrt( hc[dh] * ho[dh] );
			disCB += sqrt( hc[dh] * hb[dh] );
		}
		dis = (1 - disCO) - (1 - disCB);
		return	dis;
	}

};


/************************************************************************/
/* 目标定义                                                             */
/************************************************************************/
class	Obj
{
public:
	Obj								( const CvRect& rect , IplImage* pImg ){
		// 设置好框的位置和显示的颜色
		_rect		= rect;
		_color[0]	= 255;
		_color[1]	= 0;
		_color[2]	= 0;
		// 分配 1 当前和粒子的直方图(只是其中一个粒子的直方图),
		//		2 粒子的权重
		//		3 粒子的位置和尺度 ( 当前 )
		//		4 粒子的位置和尺度 ( 过去 )
		_curHist		= new( double[Constant::_Resolution] );
		_curHistU		= new( double[Constant::_Resolution] );
		_curBGHist		= new( double[Constant::_Resolution] );
		_particleHist	= new( double[Constant::_Resolution] );
		_weight			= new( double[Constant::_particleNum] );
		_tmpWeight		= new( double[Constant::_particleNum] );
		_cumWeightSum	= new( double[Constant::_particleNum] );
		_ind			= new( int[Constant::_particleNum] );
		_curParticlePS	= new( double[Constant::_particleNum*3] );
		_pastParticlePS	= new( double[Constant::_particleNum*3] );
		_predParticlePS	= new( double[Constant::_particleNum*3] );
		_noise			= new( double[Constant::_particleNum*3] );
		_tmpCurParticlePS	= new( double[Constant::_particleNum*3] );
		_tmpPastParticlePS	= new( double[Constant::_particleNum*3] );

		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			_weight[dp] = _tmpWeight[dp] = 1.0 / Constant::_particleNum;
		}

		_hsvLambda		= -20.0;
		_updateRatioHist= 0.05;

		IplImage*	pHSV = cvCreateImage( cvSize(pImg->width,pImg->height) , pImg->depth , pImg->nChannels );
		cvCvtColor( pImg , pHSV , CV_RGB2HSV );
		CalHist::CalHSVHist( pHSV , _rect , _curHist );
		cvReleaseImage( &pHSV );

	}
	~Obj							( ){
		if( NULL != _curHist ){
			delete _curHist;
			_curHist = NULL;
		}
		if( NULL != _curHistU ){
			delete _curHistU;
			_curHistU = NULL;
		}
		if( NULL != _curBGHist ){
			delete _curBGHist;
			_curBGHist = NULL;
		}
		if( NULL != _particleHist ){
			delete _particleHist;
			_particleHist = NULL;
		}
		if( NULL != _weight ){
			delete _weight;
			_weight = NULL;
		}
		if( NULL != _tmpWeight ){
			delete _tmpWeight;
			_tmpWeight = NULL;
		}
		if( NULL != _cumWeightSum ){
			delete _cumWeightSum;
			_cumWeightSum = NULL;
		}
		if( NULL != _ind ){
			delete _ind;
			_ind = NULL;
		}
		if( NULL != _curParticlePS ){
			delete _curParticlePS;
			_curParticlePS = NULL;
		}
		if( NULL != _pastParticlePS ){
			delete _pastParticlePS;
			_pastParticlePS = NULL;
		}
		if( NULL != _tmpCurParticlePS ){
			delete _tmpCurParticlePS;
			_tmpCurParticlePS = NULL;
		}
		if( NULL != _tmpPastParticlePS ){
			delete _tmpPastParticlePS;
			_tmpPastParticlePS = NULL;
		}
		if( NULL != _predParticlePS ){
			delete _predParticlePS;
			_predParticlePS = NULL;
		}
		if( NULL != _noise ){
			delete _noise;
			_noise = NULL;
		}
	}
public:
	// 初始化
	void							Initialize		( void ){
		CvPoint	center = GlobalMethod::GblGetCenter( _rect );
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			_curParticlePS[dp*3+0]	= gRand.GaussRand()*0.5 + center.x;
			_curParticlePS[dp*3+1]	= gRand.GaussRand()*0.5 + center.y;
			_curParticlePS[dp*3+2]	= 1.0;
			_pastParticlePS[dp*3+0] = _curParticlePS[dp*3+0];
			_pastParticlePS[dp*3+1] = _curParticlePS[dp*3+1];
			_pastParticlePS[dp*3+2] = 1.0;
		}
	}

	// 生成位置信息
	void							GenParticle		( void ){
		CvPoint	center = GlobalMethod::GblGetCenter( _rect );
		// 得到噪声
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			_noise[dp*3+0] = gRand.GaussRand()*Constant::_priorNoiseWidth;
			_noise[dp*3+1] = gRand.GaussRand()*Constant::_priorNoiseHeight;
			_noise[dp*3+2] = gRand.GeneralRandP5()*Constant::_priorNoiseScale;
		}
		// 选择下一个位置
		double	totalPrior = 0.0;
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			if( gRand.GeneralRand( ) > 0.5 ){
				for( int ind=0; ind<3; ++ind ){
					_predParticlePS[dp*3+ind] = _curParticlePS[dp*3+ind] + _noise[dp*3+ind];
				}
			} else {
				for( int ind=0; ind<3; ++ind ){
					_predParticlePS[dp*3+ind] = 2*_curParticlePS[dp*3+ind] - _pastParticlePS[dp*3+ind] + _noise[dp*3+ind];
				}
			}
			_predParticlePS[dp*3+2] = 1.0 + _noise[dp*3+2] ;
			// 尺度的范围在0.95～1.05之间
			_predParticlePS[dp*3+2] = max( 0.95 , _predParticlePS[dp*3+2] );
			_predParticlePS[dp*3+2] = min( 1.05 , _predParticlePS[dp*3+2] );
		}
		// 更新位置信息
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			for( int ind=0; ind<3; ++ind ){
				_pastParticlePS[dp*3+ind] = _curParticlePS[dp*3+ind];
				_curParticlePS[dp*3+ind]  = _predParticlePS[dp*3+ind];
			}
		}
	}

	// 先验概率的生成
	void							GenObsvLkhd		( const IplImage* pImg ){
		CvRect	searchRect;
		double	minWeight = 100.0 , totalWeight = 0;
		// 计算采样点和当前直方图之间的距离
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			searchRect = GlobalMethod::GblLocalRect(  pImg 
													, _rect 
													, _curParticlePS[dp*3+0] 
													, _curParticlePS[dp*3+1] 
													, _curParticlePS[dp*3+2] );
			CalHist::CalHSVHist( pImg , searchRect , _particleHist );
			_weight[dp] = CalHist::CalHistCoffDist( _curHist , _particleHist );
			if( _weight[dp] < minWeight ){
				minWeight = _weight[dp];
			}
		}
		// 计算权重 weight = exp( -lambda*(dis-minDis) )
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			_weight[dp] = exp( _hsvLambda*(_weight[dp]-minWeight) );
			totalWeight+= _weight[dp];
		}
		// 归一化权重
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			_weight[dp] = _weight[dp] / totalWeight;
		}
	}

	// 更新重要性
	void							UpdateImpWeight	( const IplImage* pImg ){
		// 计算出积分权重
		_cumWeightSum[0] = _weight[0];
		for( int dp=1; dp<Constant::_particleNum; ++dp ){	
			_cumWeightSum[dp] = _cumWeightSum[dp-1] + _weight[dp];
		}
		// 计算粒子位置
		for( int dp=1; dp<Constant::_particleNum; ++dp ){
			_ind[dp] = 0;
		}
		int		ind = 0;
		double	u	= gRand.GeneralRand()/Constant::_particleNum , uInd;
		for( int dp=0; dp<Constant::_particleNum; ++dp ){	
			uInd = u + double(dp)/Constant::_particleNum;
			while( uInd > _cumWeightSum[ind] ){
				ind = ind + 1;
			}
			_ind[dp] = ind;
		}
		
		// 重新分配粒子
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			for( int ind=0; ind<3; ++ind ){
				_tmpCurParticlePS[dp*3+ind]  = _curParticlePS[_ind[dp]*3+ind];
				_tmpPastParticlePS[dp*3+ind] = _pastParticlePS[_ind[dp]*3+ind];
			}
			_tmpWeight[dp] = _weight[_ind[dp]];
		}
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			for( int ind=0; ind<3; ++ind ){
				_curParticlePS[dp*3+ind]  = _tmpCurParticlePS[dp*3+ind];
				_pastParticlePS[dp*3+ind] = _tmpPastParticlePS[dp*3+ind];
			}
			_weight[dp] = _tmpWeight[dp];
		}
		// 计算位置和尺度
		double	x = 0.0 , y = 0.0 , w = 0.0 , h = 0.0;
		for( int dp=0; dp<Constant::_particleNum; ++dp ){
			x += _curParticlePS[dp*3+0] - _rect.width / 2.0;
			y += _curParticlePS[dp*3+1] - _rect.height/ 2.0;
			w += _curParticlePS[dp*3+2] * _rect.width;
			h += _curParticlePS[dp*3+2] * _rect.height;
		}
		x = x / Constant::_particleNum;		y = y / Constant::_particleNum;
		w = w / Constant::_particleNum;		h = h / Constant::_particleNum;
		_rect.x		= cvFloor( x+0.5 );		_rect.y		= cvFloor( y+0.5 );
		_rect.width = cvFloor( w+0.5 );		_rect.height= cvFloor( h+0.5 );
		GlobalMethod::GblRectCheckBound( &_rect , pImg->width , pImg->height );
	}

	// 模板更新
	void							UpdateTemplate( IplImage* pImg ){
		CalHist::CalHSVHist( pImg , _rect , _curHistU );
		for( int dh=0; dh<Constant::_Resolution; ++dh ) {
			_curHist[dh] = _curHist[dh]*(1-_updateRatioHist)+_curHistU[dh]*_updateRatioHist;
		}
	}

public:
	CvRect							_rect;				// 跟踪框
	int								_color[3];			// 跟踪框的颜色（用于画）
	double*							_curHist;			// 模板目标直方图
	double*							_curHistU;			// 模板目标直方图(用于更新)
	double*							_curBGHist;			// 模板背景直方图
	double*							_particleHist;		// 每个粒子的直方图
	double*							_weight;			// 粒子权重
	double*							_tmpWeight;			// 粒子权重(用于重采样)
	double*							_cumWeightSum;		// 粒子权重积分(用于重采样)
	int*							_ind;				// 粒子采样Index
	double*							_curParticlePS;		// 当前粒子状态
	double*							_pastParticlePS;	// 过去粒子状态
	double*							_tmpCurParticlePS;
	double*							_tmpPastParticlePS;
	double*							_predParticlePS;	// 预测粒子的状态
	double*							_noise;				// 噪声
	double							_hsvLambda;			// lambda值（-20）
	double							_updateRatioHist;	// 直方图的更新率
};

/************************************************************************/
/* 目标列表的定义                                                       */
/************************************************************************/
class	TrackObjList
{
public:
	TrackObjList					( IplImage* pImg )
		:	_version( "1.0.0.0" )
	{
		// 读取图像
		if( NULL != pImg ){
			_imWidth = pImg->width;
			_imHeight= pImg->height;
		} else {
			exit( 0 );
		}
		// 物体列表清空
		_objList.clear( );
		// 输出提示信息
		cout	<< "目前的版本是:\t\t"		<< _version		<< endl;
		cout	<< "图像的宽度是:\t\t"		<< _imWidth		<< endl;
		cout	<< "图像的高度是:\t\t"		<< _imHeight	<< endl;
		cout	<< "粒子的数目是:\t\t"		<< Constant::_particleNum	<< endl;
		// 帧号设置为0
		_frameInd	= 0;
	}

	~TrackObjList					( ){
		for( vector<Obj*>::iterator	iter=_objList.begin(); iter<_objList.end(); ++iter ){
			delete	(*iter);
		}
		_objList.clear( );
	}

public:
	// 检测新物体
	void							DetectObj		( IplImage* pImg ){
		++ _frameInd;
		if( _frameInd == 1 ){
			// 初始化物体的位置
			CvRect	rect = cvRect( 147 , 92 , 20  , 30 );
			Obj*	pObj = new Obj( rect , pImg );
			pObj->Initialize( );
			_objList.push_back( pObj );
		}
	}

	// 粒子滤波过程( 也就是更新过程 )
	void							Update			( IplImage* pImg ){
		IplImage*	pHSV = cvCreateImage( cvSize(pImg->width,pImg->height) , pImg->depth , pImg->nChannels );
		cvCvtColor( pImg , pHSV , CV_RGB2HSV );
		for( vector<Obj*>::iterator	iter=_objList.begin(); iter<_objList.end(); ++iter ){
			(*iter)->GenParticle( );
			(*iter)->GenObsvLkhd( pHSV );
			(*iter)->UpdateImpWeight( pHSV );
			(*iter)->UpdateTemplate( pHSV );
		}
		cvReleaseImage( &pHSV );
	}
	// 显示结果
	void							DrawImg			( IplImage* pImg ){
		for( vector<Obj*>::iterator	iter=_objList.begin(); iter<_objList.end(); ++iter ){
			cvRectangle(  pImg 
						, cvPoint( (*iter)->_rect.x , (*iter)->_rect.y ) 
						, cvPoint( (*iter)->_rect.x + (*iter)->_rect.width , (*iter)->_rect.y + (*iter)->_rect.height )
						, CV_RGB((*iter)->_color[0],(*iter)->_color[1],(*iter)->_color[2]) 
						, 2
						, CV_AA 
						, 0 );
		}
	}
public:
	vector<Obj*>					_objList;					// 物体列表
	int								_imWidth;					// 图像宽度
	int								_imHeight;					// 图像高度
	string							_version;					// 版本号
	int								_frameInd;					// 帧号
};

/************************************************************************/
/* 图片序列读取                                                         */
/************************************************************************/
class ImListReader 
{
public:
	ImListReader( )
		:	m_StartNumber( 1 )
		,	m_EndNumber( 0 )
	{
	}

	~ImListReader( )
	{
		if( NULL != m_pImg ){
			cvReleaseImage( &m_pImg );
			m_pImg = NULL;
		}
	}

	bool			OpenVideo( string filePathHead , string ctrlFileName ) {
		string filePath = filePathHead;
		filePath.append( ctrlFileName );
		ifstream inputStream;
		inputStream.open( filePath.c_str() , ios::in );
		if( inputStream.fail() ) {
			return false;
		} else {
			/***
			 * 读取控制文件信息
			 ***/
			string line;
			getline( inputStream , line , '\n' );
			m_StartNumber = stoi( line );
			m_CurNumber = m_StartNumber;
			getline( inputStream , line , '\n' );
			m_EndNumber = stoi( line );
			getline( inputStream , line , '\n' );
			m_FileFormat = filePathHead;
			m_FileFormat.append( line );
			/***
			 * 关闭文件
			 ***/
			inputStream.close( );
			return true;
		}
	}

	bool			ReadFrame( ){

		if ( m_CurNumber > m_EndNumber ) {
			return false;
		}

		string imgFilePath = format( m_FileFormat.c_str() , m_CurNumber );
		m_pImg = cvLoadImage( imgFilePath.c_str( ) );
		if ( NULL == m_pImg ) {
			return false;
		}

		if ( 1 == m_pImg->nChannels ){
			return false;
		}

		++ m_CurNumber;
		return true;
	}

public:
	IplImage*		m_pImg;					// 当前图像(必须为3通道)

protected:
	int				m_StartNumber;			// 初始帧号
	int				m_EndNumber;			// 末尾帧号
	int				m_CurNumber;			// 当前帧号
	string			m_FileFormat;			// 文件读取文件路径

private:
	std::string format(const char *fmt, ... ) 
	{ 
		std::string strResult="";
		if (NULL != fmt) {
			va_list marker = NULL;            
			va_start(marker, fmt);                           // 初始化变量参数 
			size_t nLength = _vscprintf(fmt, marker) + 1;    // 获取格式化字符串长度
			std::vector<char> vBuffer(nLength, '\0');        // 创建用于存储格式化字符串的字符数组
			int nWritten = _vsnprintf_s(&vBuffer[0], vBuffer.size(), nLength, fmt, marker);
			if (nWritten>0) {
				strResult = &vBuffer[0];
			}            
			va_end(marker);                                    // 重置变量参数
		}
		return strResult; 
	}
};



#endif