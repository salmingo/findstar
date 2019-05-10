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

}

// 判定参考星是否位于搜索范围内
bool ACatalog::inner_zone(double ra, double dec) {
	return false;
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

	return false;
}

//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */
