/*!
 * class ATimeSpace 常用天文时空相关转换接口
 * Version：0.2
 *
 * 功能列表:
 * 1. 赤经/赤纬格式转换: 字符串<->实数
 * 2. 时间/角度格式转换: 实数<->时/度, 分, 秒
 * 3. 儒略日和修正儒略日
 * 4. 平恒星时和真恒星时
 *
 * 创建时间: 2012年3月1日
 * 完成时间: 2012年3月10日
 * 作者: 卢晓猛, lxm@nao.cas.cn
 *
 * 暂缺项（2015-02-15）
 * <1> 力学时和世界协调时的相互转换
 * <2> 天体的升起、过中天及降落时间
 * <3> 大气折射
 */

#ifndef _ATIMESPACE_H_
#define _ATIMESPACE_H_

namespace AstroUtil {
///////////////////////////////////////////////////////////////////////////////
// ATimeSpace输出数据内部索引, 用于减少重复计算以提高速度
#define DAT			0		//< 国际原子时和世界调协时之差=TAI - UTC
#define	JD			1		//< 儒略日
#define MJD			2		//< 修正儒略日
#define EPOCH		3		//< 历元
#define JC			4		//< 儒略世纪
#define GMST		5		//< 格林尼治0时的平恒星时
#define GAST		6		//< 格林尼治0时的真恒星时
#define LMST		7		//< 测站位置的平恒星时
#define LAST		8		//< 测站位置的真恒星时
#define NL			9		//< 章动的经度项
#define NO			10		//< 章动的纬度项
#define MLA			11		//< 月球纬度参数, MoonLatArg
#define MLAN		12		//< 月球升交点经度, MoonLongAscendingNode
#define MME			13		//< 月球到太阳的平角距离, MoonMeanElongation
#define MA			14		//< 地球绕日公转的平近点角, MeanAnomaly
#define MMA			15		//< 月球的平近点角, MoonMeanAnomaly
#define MEO			16		//< 平黄赤交角, MeanEclipObliquity
#define TEO			17		//< 真黄赤交角, TrueEclipObliquity

#define LEN_ATS		30		//< ATimeSpace数据缓冲区长度

class ATimeSpace
{
public:
	ATimeSpace();
	virtual ~ATimeSpace();

private:
	double	m_val[LEN_ATS];		//< 为避免重复计算, 类函数计算结果的临时缓存区
	bool	m_stat[LEN_ATS];	//< 为避免重复计算, 类函数计算结果的更新状态
	double	m_lgt;			//< 地理经度, 量纲: 弧度
	double	m_lat;			//< 地理纬度, 量纲: 弧度
	double	m_alt;			//< 海拔高度, 量纲: 米

public:
	/*!
	 * @brief 设置时间, 以后的计算都基于该时间执行
	 * @param[in] year   年
	 * @param[in] month  月
	 * @param[in] day    日, 在月中的日数, 从1开始, 即日历中的日期
	 * @param[in] hour   时, 当天0时为零点的小时数
	 */
	void SetTime(int year, int month, int day, double hour);
	/*!
	 * @brief 地理位置, 以后的相关计算都基于该位置执行
	 * @param[in] lgt 地理经度, 量纲: 角度
	 * @param[in] lat 地理纬度, 量纲: 角度
	 * @param[in] alt 海拔高度, 量纲: 米
	 */
	void SetSite(double lgt, double lat, double alt);

// Function
public:
	/*!
	 * \brief Transfer R.A. from string to double
	 * \param[in] pszVal R.A. in string style
	 * \param[out] val    R.A. in double, in hour
	 * \return
	 * if pszVal format is right then transfer it to be double and return true, or else return false
	 * \note
	 * pszVal could be the following style:
	 * HH:MM:SS.SS, colon could be omitted or replaced by space, and both MM and SS could be omitted too.
	 */
	bool RAStr2Dbl(const char *pszVal, double &val);
	/*!
	 * \brief Transfer DEC. from string to double
	 * \param[in] pszVal DEC. in string style
	 * \param[out] val    DEC. in double, in degree
	 * \return
	 * if pszVal format is right then transfer it to be double and return true, or else return false
	 * \note
	 * pszVal could be the following style:
	 * sDD:MM:SS.SS, colon could be omitted or replaced by space, and both MM and SS could be omitted too.
	 * 's' means plus or minus sign
	 */
	bool DECStr2Dbl(const char *pszVal, double &val);
	/*!
	 * \brief Transfer R.A. from double to string
	 * \param[in]    val      R.A. in double, in hour
	 * \param[out] pszVal R.A. in string style
	 * \return
	 * regulate val to be between 0 and 24 and transfer it to be string style
	 * \note
	 * pszVal could be the following style:
	 * HH:MM:SS.SS, colon could be omitted or replaced by space, and both MM and SS could be omitted too.
	 */
	void RADbl2Str(const double val, char *pszVal);
	/*!
	 * \brief Transfer DEC. from double to string
	 * \param[in]    val      DEC. in double, in hour
	 * \param[out] pszVal DEC in string style
	 * \return
	 * regulate val to be between -90 and +90 and transfer it to be string style
	 * \note
	 * pszVal could be the following style:
	 * sDD:MM:SS.SS, colon could be omitted or replaced by space, and both MM and SS could be omitted too.
	 * 's' means plus or minus sign
	 */
	void DECDbl2Str(const double val, char *pszVal);
	/*!
	 * \brief resolve double hour to hour, minute and second
	 * \param[in] hour double time, in hour
	 * \param[out] hh  hour
	 * \param[out] mm  minute
	 * \param[out] ss  second
	 */
	void Hour2HMS(const double hour, int &hh, int &mm, double &ss);
	/*!
	 * \brief construct hour, minute and second to double hour
	 * \param[in] hh  hour
	 * \param[in] mm  minute
	 * \param[in] ss  second
	 * \return
	 * double hour
	 */
	double HMS2Hour(const int hh, const int mm, const double ss);
	/*!
	 * \brief resolve double degree to degree, minute and second
	 * \param[in] deg double degree, in degree
	 * \param[out] dd  degree
	 * \param[out] mm  minute
	 * \param[out] ss  second
	 * \param[out] sign plus or minus sign. if deg is less than 0 then sign is -1, or else 1
	 */
	void Deg2DMS(const double deg, int &dd, int &mm, double &ss, int &sign);
	/*!
	 * \brief construct degree, minute and second to double degree
	 * \param[in] sign plus or minus sign. if deg is less than 0 then sign is -1, or else 1
	 * \param[in] dd  degree
	 * \param[in] mm  minute
	 * \param[in] ss  second
	 * \return
	 * double degree
	 */
	double DMS2Deg(const int sign, const int dd, const int mm, const double ss);

public:
	/*!
     * @brief whether this is a leap year
	 */
	bool LeapYear(int year);
	/*!
	 * @brief date is the nth day of this year
	 */
	int DayOfYear(int year, int month, int day);
	/*!
	 * @brief the date corresponding to the nth day of this year
	 */
	void Day2Date(int year, int nth, int &month, int &day);
	/*!
	 * @brief the day of a week
	 */
	int DayOfWeek(int year, int month, int day);

public:
	/*!
	 * @brief 计算输入时间对应的历元
	 * @param[in] mjd 修正儒略日
	 * @return
	 * 设置时间或儒略日对应的历元
	 */
	double Epoch();
	double Epoch(double mjd);
    /*!
	 *  @brief 根据历元计算对应的儒略世纪
	 *  @param[in] ep 历元
	 *  @return
	 *  输入时间对应的儒略世纪. 从J2000开始计算
	 */
	double Epoch2JC();
	double Epoch2JC(double ep);
	/*!
	 * @brief 由历书时计算儒略日
	 *  \param[in] year  年
	 *  \param[in] month 月
	 *  \param[in] day   日
	 *  \param[in] hour  当日小时数
	 *  \param[in] eph   历书标记. 1: 格里高里历; 2: 儒略历
	 * @return
	 * 公元历对应的儒略日
	 */
	double JulianDay(int year, int month, int day, double hour, int eph = 1);
    /*!
	 *  \brief 计算修正儒略日
	 *  \param[in] year 年
	 *  \param[in] month 月
	 *  \param[in] day 日
	 *  \param[in] hour 时
	 *  \param[in] jd   儒略日
	 *  \return 年月日时对应的修正儒略日
	 **/
	double ModifiedJulianDay(int year, int month, int day, double hour);
	double ModifiedJulianDay(double jd);
	/*!
	 * @brief 由儒略日计算对应的公元历
	 * @param jd    儒略日
	 * @param year  年
	 * @param month 月
	 * @param day   日
	 * @param hour  当日小时数
	 */
	void JD2YMDH(double jd, int &year, int &month, int &day, double &hour);
	/*!
	 * @brief 时区力学时转换为国际原子时
	 * @param[in] tt 力学时, 量纲: 天
	 * @return
	 * 国际原子时, 量纲: 天
	 */
	double TT2TAI(double tt);
	/*!
	 * @brief 国际原子时转换为时区力学时
	 * @param[in] tai 原子时, 量纲: 天
	 * @return
	 * 力学时, 量纲: 天
	 */
	double TAI2TT(double tai);
	/*!
	 * @brief 计算在某一世界协调时, 国际原子时和世界协调时的差异
	 * @param[in]  year   年
	 * @param[in]  month  月
	 * @param[in]  day    日
	 * @param[in]  hour   当日小时数
	 * @param[out] dat    闰秒, 量纲: 天
	 * @return
	 *   1 -- 离有效期较远, 值得怀疑其准确性
	 *   0 -- 结果正确
	 *  -1 -- 年, 错误
	 *  -2 -- 月, 错误
	 *  -3 -- 日, 错误
	 *  -4 -- 时, 错误
	 */
	int DeltaAT(int year, int month, int day, double hour, double &dat);

public:
	/*!
	 * @brief 根据修正儒略日计算格林威治0时平恒星时.
	 * @param[in] mjd 修正儒略日. 修正儒略日以J2000为零点
	 * @return
	 * 以弧度为单位的平恒星时
	 */
	double GMST0();
	double GMST0(double mjd);
	/*!
	 * @brief 根据修正儒略日计算格林威治0时真恒星时
	 * @param[in] mjd 修正儒略日. 修正儒略日以J2000为零点
	 * @return
	 * 真恒星时, 量纲: 弧度
	 */
	double GAST0();
	double GAST0(double mjd);
	/*!
	 * @brief 根据修正儒略日计算本地平恒星时.
	 * @param[in] mjd 修正儒略日. 修正儒略日以J2000为零点
	 * @return 以弧度为单位的本地平恒星时
	 */
	double LocalMeanSiderialTime();
	double LocalMeanSiderialTime(double mjd);
	/*!
	 * @brief 根据修正儒略日计算本地真恒星时
	 * @param[in] mjd 修正儒略日. 修正儒略日以J2000为零点
	 * @return
	 * 本地真恒星时, 量纲: 弧度
	 */
	double LocalSiderialTime();
	double LocalSiderialTime(double mjd);

	/*!
	 * @brief 自行改正. 计算自历元1至历元2的自行改正. 历元2是SetTime()对应的数值
	 * @param[in]  ep1    历元1
	 * @param[in]  ra0    历元1对应的平位置赤经, 量纲: 弧度
	 * @param[in]  dec0   历元1对应的平位置赤纬, 量纲: 弧度
	 * @param[in]  pmra   赤经自行, 量纲: 毫角秒/年
	 * @param[in]  pmdec  赤纬自行, 量纲: 毫角秒/年
	 * @param[out] ra     改正后赤经, 量纲: 弧度
	 * @param[out] dec    改正后赤纬, 量纲: 弧度
	 */
	void CorrectPM(double ep1, double ra0, double dec0, double pmra, double pmdec, double &ra, double &dec);

	/*!
	 * @brief 历元转换, 将历元1对应的赤道平位置转换为历元2对应的赤道平位置. 历元2是SetTime()对应的数值
	 * @param[in]  ep1  历元1
	 * @param[in]  ra1  历元1对应的平位置赤经, 量纲: 弧度
	 * @param[in]  dec1 历元1对应的平位置赤纬, 量纲: 弧度
	 * @param[out] ra2  历元2对应的平位置赤经, 量纲: 弧度
	 * @param[out] dec2 历元2对应的平位置赤纬, 量纲: 弧度
	 * @note
	 * 历元转换之前, (ra2, dec1)应先做自行改正
	 */
	void EpochTransform(double ep1, double ra1, double dec1, double &ra2, double &dec2);

public:
	/*!
	 * @brief 月球纬度参数
	 * @param[in] ep 历元
	 * @return
	 * 量纲: 弧度
	 */
	double MoonLatArg();
	double MoonLatArg(double ep);
	/*!
	 * @brief 月球平轨道和黄道的升交点经度
	 * @param[in] ep 历元
	 * @return
	 * 升交点经度, 量纲: 弧度
	 */
	double MoonLongAscendingNode();
	double MoonLongAscendingNode(double ep);
	/*!
	 * @brief 章动项的经度分量
	 * @param[in] ep 历元
	 * @return
	 * 章动项经度分量, 量纲: 弧度
	 */
	double NutationLongitude();
	double NutationLongitude(double ep);
	/*!
	 * @brief 章动项的纬度分量
	 * @param[in] ep 历元
	 * @return
	 * 章动项纬度分量, 量纲: 弧度
	 */
	double NutationObliquity();
	double NutationObliquity(double ep);
	/*!
	 * @brief 计算平黄赤交角
	 * @param[in] ep 历元
	 * @return
	 * 黄道和赤道交角的平均值, 量纲: 弧度
	 * @note
	 * I7 3GHz CPU, 计算时间约: 0.33微秒
	 */
	double MeanEclipObliquity();
	double MeanEclipObliquity(double ep);
	/*!
	 * @brief 计算真黄赤交角
	 * @param[in] ep 历元
	 * @return
	 * 黄道和赤道交角的真值, 平值加章动项, 量纲: 弧度
	 * @note
	 * I7 3GHz CPU, 计算时间约: 1.03微秒
	 */
	double TrueEclipObliquity();
	double TrueEclipObliquity(double ep);
	/*!
	 * @brief 月球到太阳的平均角距离, 等效于月相
	 * @param[in] ep 历元
	 * @return
	 * 月球到太阳的平均角距离, 量纲: 角度.
	 * @note
	 * 0对应新月, 180对应满月, 90和270分别对应上弦月和下弦月
	 */
	double MoonMeanElongation();
	double MoonMeanElongation(double ep);
	/*!
	 * @brief 地球绕日公转的平近点角
	 * @param[in] ep 历元
	 * @return
	 * 平近点角, 量纲: 弧度
	 * @note
	 * 平近点角: 轨道上的物体在辅助圆上相对于中心点的角度. 参考"平近点角定义"
	 */
	double MeanAnomaly();
	double MeanAnomaly(double ep);
	/*!
	 * @brief 月球绕地公转的平近点角
	 * @param[in] ep 历元
	 * @return
	 * 平近点角, 量纲: 弧度
	 */
	double MoonMeanAnomaly();
	double MoonMeanAnomaly(double ep);

public:
	/*!
	 * @brief 黄道坐标转换到赤道坐标
	 * @param[in]  l   黄经, 量纲: 弧度
	 * @param[in]  b   黄纬, 量纲: 弧度
	 * @param[out] ra  赤经, 量纲: 弧度
	 * @param[out] dec 赤纬, 量纲: 弧度
	 * @note
	 * 黄赤转换需要当前时间
	 */
	void Eclip2Eq(double l, double b, double &ra, double &dec);
	/*!
	 * @brief 赤道坐标转换到黄道坐标
	 * @param[in]  ra  赤经, 量纲: 弧度
	 * @param[in]  dec 赤纬, 量纲: 弧度
	 * @param[out] l   黄经, 量纲: 弧度
	 * @param[out] b   黄纬, 量纲: 弧度
	 * @note
	 * 黄赤转换需要当前时间
	 */
	void Eq2Eclip(double ra, double dec, double &l, double &b);
	/*!
	 * @brief 地平坐标转换到赤道坐标
	 * @param[in]  azi 方位角, 量纲: 弧度
	 * @param[in]  alt 俯仰角, 量纲: 弧度
	 * @param[out] ha  时角, 量纲: 弧度
	 * @param[out] dec 赤纬, 量纲: 弧度
	 * @note
	 * 地平\赤道转换需要测站地理位置
	 * 方位角为南零点
	 */
	void AltAzi2Eq(double azi, double alt, double &ha, double &dec);
	/*!
	 * @brief 赤道坐标转换到地平坐标
	 * @param[in]  ha  时角, 量纲: 弧度
	 * @param[in]  dec 赤纬, 量纲: 弧度
	 * @param[out] azi 方位角, 量纲: 弧度
	 * @param[out] alt 俯仰角, 量纲: 弧度
	 * @note
	 * 地平\赤道转换需要测站地理位置
	 * 方位角为南零点
	 */
	void Eq2AltAzi(double ha, double dec, double &azi, double &alt);

public:
	/**
	 * \brief 计算地平式望远镜抵消地平自转, 消转器的转动弧度
	 * \param[in] ha   时角, 量纲: 弧度
	 * \param[in] dec  赤纬, 量纲: 弧度
	 * \param[in] mode 卡塞格林焦点, mode='N'
	 *                 耐氏焦点, 东边焦点mode='+'
	 *                 西边焦点mode='-'
	 * \return 抵消地球自转需要的弧度
	 **/
	double ParAngle(double ha, double dec, char mode);

};
///////////////////////////////////////////////////////////////////////////////
}

#endif /* _ATIMESPACE_H_ */
