/*
 * @file ACatalog.cpp 天文星表访问共用接口
 */

#include "ADefine.h"
#include "ACatalog.h"

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
ACatalog::ACatalog() {
	ra_      = dec_     = r_ = 0.0;
	ra_min_  = ra_max_  = 0.0;
	dec_min_ = dec_max_ = 0.0;
}

ACatalog::~ACatalog() {

}

// 初始化搜索参数
void ACatalog::init_param() {
//	double sf = sin(r_);
//	double cd = cos(dec_);
//	double d;

//	if ((spdmin = (int) ((dec_ - r_)) < 0) spdmin = 0;
//			if ((spdmax = (int) ((dec_ + r_ + 90.) * MILLISEC)) > MILLISEC180) spdmax = MILLISEC180 - 1;
//	if (sf < cd) {
//		d = asin(sf / cd) * RtoG;
//		if ((ramin = (int) ((ra - d) * MILLISEC)) < 0)            ramin += MILLISEC360;
//		if ((ramax = (int) ((ra + d) * MILLISEC)) >= MILLISEC360) ramax -= MILLISEC360;
//	}
//	else {
//		ramin = 0;
//		ramax = MILLISEC360 - 1;
//	}
//	if (ramin > ramax) ramax += MILLISEC360;
}

// 判定参考星是否位于搜索范围内
bool ACatalog::inner_zone(double ra, double dec) {
	double v = cos(dec) * cos(dec_) * cos(ra - ra_) + sin(dec) * sin(dec_);
	return acos(v) <= r_;
}

// 设置星表根路径
void ACatalog::SetPath(const char *path) {
	pathCat_ = path;
}

// 从星表中查找视场范围内的参考星
bool ACatalog::FindStar(double ra, double dec, double r) {
	ra_  = ra * D2R;
	dec_ = dec * D2R;
	r_   = r * D2R;
	init_param();

	return true;
}

//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */
