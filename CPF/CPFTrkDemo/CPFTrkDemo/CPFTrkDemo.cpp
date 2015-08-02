// CPFTrkDemo.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include <opencv2\opencv.hpp>  
#include <iostream>  
#include <string>  
using namespace cv;  
using namespace std;  

#include "CPF.hpp"

int _tmain(int argc, _TCHAR* argv[])
{

	/************************************************************************/
	/* 视频文件控制信息                                                     */
	/************************************************************************/
	string filePathHead = "./Data/";
	string ctrlFileName = "Video.txt";

	/************************************************************************/
	/* 读取视频文件的控制信息（打开文件）                                   */
	/************************************************************************/
	ImListReader imListReader;
	if ( false == imListReader.OpenVideo( filePathHead , ctrlFileName ) ) {
		return 0;
	}
	/************************************************************************/
	/* 初始化                                                               */
	/************************************************************************/
	if ( false == imListReader.ReadFrame() ){
		return 0;
	}
	TrackObjList	trackObjList( imListReader.m_pImg );

	while( true == imListReader.ReadFrame() ){

		cout << "正在处理第" << trackObjList._frameInd << "帧!" << endl;
		trackObjList.DetectObj( imListReader.m_pImg );
		trackObjList.Update( imListReader.m_pImg );
		trackObjList.DrawImg( imListReader.m_pImg );
		cvShowImage( "CPF Demo", imListReader.m_pImg );
		cvWaitKey( 1 );

	}

	return 0;
}

