/**
 Name        : ADefine.h 常量及宏定义声明文件
 Author      : Xiaomeng Lu
 Version     : 0.1
 Date        : Oct 13, 2012
 Last Date   : Oct 13, 2012
 Description : 天文数字图像处理中常用的常量及宏定义
 **/
#ifndef _ADEFINE_H_
#define _ADEFINE_H_

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <vector>
#include <string>
#include <algorithm>
#include <typeinfo>

#if !defined(WIN32) && !defined(_WINDOWS_)	// Linux平台下使用该类库
#include <unistd.h>
#else
#include <io.h>								// Windows平台使用该类库
#endif

using namespace std;

namespace AstroUtil
{
///////////////////////////////////////////////////////////////////////////////
/********************************** 宏定义 **********************************/
#if !defined(BZERO)	// 内存清零
    #define BZERO(ptr, szbyte)	memset(ptr, 0, szbyte)
#endif

#if !defined(BYTE)
    typedef unsigned char BYTE;
#endif

#if !defined(USHORT)
    typedef unsigned short USHORT;
#endif

#if !defined(UINT)
    typedef unsigned int UINT;
#endif

#ifndef WORD
    typedef unsigned short WORD;
#endif

// 圆周率相关常量
#if !defined(API)	// 圆周率
#define API		3.141592653589793
#endif
#define PI45	0.785398163397448
#define PI90	1.570796326794897
#define PI180	API
#define PI270	4.712388980384690
#define PI360	6.283185307179586
#define GtoR	0.0174532925199433		// 使用乘法, 角度转换为弧度的系数
#define RtoG	57.295779513082323		// 使用乘法, 弧度转换为角度的系数
#define StoR	4.848136811095360e-06	// 使用乘法, 角秒转换为弧度的系数
#define RtoS	2.062648062470964e+05	// 使用乘法, 弧度转换为角秒的系数
#define ASEC	StoR			// 角秒对应的弧度
// 时间相关常量
#define MILLISEC		3600000			// 1度对应的毫角秒
#define MILLISEC5		18000000		// 5度对应的毫角秒
#define MILLISEC90		324000000		// 90度对应的毫角秒
#define MILLISEC180		648000000		// 180度对应的毫角秒
#define MILLISEC360		1296000000		// 360度对应的毫角秒

#define DAYSEC			86400			// 1天内的秒数
#define TTMTAI			32.184			// TT减TAI, 量纲: 秒

// 软件中限定最大/最小数
#define AMAX	1E30
#define AEPS 	1E-4

#define cot(x)	(tan(PI90 - (x)))	// 计算弧度对应的余切

///////////////////////////////////////////////////////////////////////////////

typedef struct param_mount_site{// 定义与测站位置有关信息
	double	lgt;	// 地理经度, 量纲: 弧度
	double	lat;	// 地理纬度, 量纲: 弧度
	double	alt;	// 海拔高度, 量纲: 米
	double	airp;	// 大气压, 量纲: 百帕
	double	temp;	// 环境温度, 量纲: 摄氏度
	int		tmc;	// 时区修正, 所使用本地时的时间与地理位置所在时区的差
}* ptr_param_mntst;

// 定义实数型'点'
typedef struct point_2d {// 两值实数点
	double x;
	double y;

public:
	point_2d() {
		x = y = 0;
	}
	point_2d(double x, double y) {
		this->x = x;
		this->y = y;
	}
	point_2d& operator=(const point_2d &pt) {
		if (this != &pt)
			memcpy(this, &pt, sizeof(point_2d));
		return *this;
	}
}* ptr_point_2d;

typedef struct point_3d {// 三值实数点
	double x;
	double y;
	double z;

public:
	point_3d() {
		x = y = z = 0;
	}
	point_3d(point_3d &pt) {
		this->x = pt.x;
		this->y = pt.y;
		this->z = pt.z;
	}
	point_3d(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	point_3d& operator=(const point_2d &pt) {
		this->x = pt.x;
		this->y = pt.y;
		this->z = 0;
		return *this;
	}
	point_3d& operator=(const point_3d &pt) {
		if (this != &pt)
			memcpy(this, &pt, sizeof(point_3d));
		return *this;
	}
}* ptr_point_3d;

typedef struct obstime {//< observation time
    int     year;       // year
    int     month;      // month
    int     day;        // day
    double  hours;      // hours in day

public:
    obstime() {
        year = month = day = 0;
        hours = -1; // hours equas to 0 represents invalid time
    }

    obstime(int year, int month, int day, double hours = 0.) {
        this->year  = year;
        this->month = month;
        this->day   = day;
        this->hours = hours;
    }

    void set(int year, int month, int day, double hours) {
        this->year  = year;
        this->month = month;
        this->day   = day;
        this->hours = hours;
    }
}* ptr_obstime;

typedef struct errlog {// error log
	string utc;		// UTC when error being triggered
	string message;	// error message

public:
	void newlog(const char *msg) {
		time_t t = time(NULL);
		utc = asctime(gmtime(&t));
		message = msg;
	}
}* ptrelog;

///////////////////////////////////////////////////////////////////////////////
};
#endif
