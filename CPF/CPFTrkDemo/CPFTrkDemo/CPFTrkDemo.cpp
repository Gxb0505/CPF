// CPFTrkDemo.cpp : �������̨Ӧ�ó������ڵ㡣
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
	/* ��Ƶ�ļ�������Ϣ                                                     */
	/************************************************************************/
	string filePathHead = "./Data/";
	string ctrlFileName = "Video.txt";

	/************************************************************************/
	/* ��ȡ��Ƶ�ļ��Ŀ�����Ϣ�����ļ���                                   */
	/************************************************************************/
	ImListReader imListReader;
	if ( false == imListReader.OpenVideo( filePathHead , ctrlFileName ) ) {
		return 0;
	}
	/************************************************************************/
	/* ��ʼ��                                                               */
	/************************************************************************/
	if ( false == imListReader.ReadFrame() ){
		return 0;
	}
	TrackObjList	trackObjList( imListReader.m_pImg );

	while( true == imListReader.ReadFrame() ){

		cout << "���ڴ����" << trackObjList._frameInd << "֡!" << endl;
		trackObjList.DetectObj( imListReader.m_pImg );
		trackObjList.Update( imListReader.m_pImg );
		trackObjList.DrawImg( imListReader.m_pImg );
		cvShowImage( "CPF Demo", imListReader.m_pImg );
		cvWaitKey( 1 );

	}

	return 0;
}

