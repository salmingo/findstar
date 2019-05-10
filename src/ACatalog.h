/*
 * @file ACatalog.h 星表访问共用接口
 * @date 2019年5月10日
 * @version 0.1
 */

#ifndef ACATALOG_H_
#define ACATALOG_H_

#include <string>
using std::string;

namespace AstroUtil {
//////////////////////////////////////////////////////////////////////////////
/*!
 * @class ACatalog 天文星表访问共用接口
 * @note
 * 按照输入参数, 从星表中查找符合条件的参考星
 */
class ACatalog {

public:
	ACatalog();
	virtual ~ACatalog();

protected:
	string pathCat_;	//< 星表路径
	/*!
	 * @member ra  视场中心赤经, 量纲: 弧度
	 * @member dec 视场中心赤纬, 量纲: 弧度
	 * @member r   视场半径, 量纲: 弧度
	 */
	double ra_;
	double dec_;
	double r_;
	/*!
	 * @member ra_min   视场覆盖区域最小赤经
	 * @member ra_max   视场覆盖区域最大赤经
	 * @member dec_min  视场覆盖区域最小赤纬
	 * @member dec_max  视场覆盖区域最大赤纬
	 */
	double ra_min_, ra_max_;	//< 赤经范围, 量纲: 角度
	double dec_min_, dec_max_;	//< 赤纬范围, 量纲: 角度

protected:
	/*!
	 * @brief 初始化搜索参数
	 */
	void init_param();
	/*!
	 * @brief 判定参考星是否位于搜索范围内
	 * @param ra   赤经, 量纲: 弧度
	 * @param dec  赤纬, 量纲: 弧度
	 * @return
	 * 判定结果. 参考星位于搜索范围内则返回true
	 */
	bool inner_zone(double ra, double dec);

public:
	/*!
	 * @brief 设置星表根路径
	 * @param path 星表数据根路径
	 */
	void SetPath(const char *path);
	/*!
	 * @brief 从星表中查找视场范围内的参考星
	 * @param ra   视场中心赤经, 量纲: 角度
	 * @param dec  视场中心赤纬, 量纲: 角度
	 * @param r    视场半径, 量纲: 角度
	 * @return
	 * 参考星查询结果
	 */
	virtual bool FindStar(double ra, double dec, double r);
};
//////////////////////////////////////////////////////////////////////////////
} /* namespace AstroUtil */

#endif /* ACATALOG_H_ */
