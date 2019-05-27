/*
 * ACatUCAC4 访问UCAC4星表封装类, 定义文件
 */

#include "ACatUCAC4.h"
#include "ATimeSpace.h"
#include "AMath.h"

namespace AstroUtil
{
///////////////////////////////////////////////////////////////////////////////
ACatUCAC4::ACatUCAC4()
	: ACatalog() {
	m_stars = NULL;
	m_asc   = NULL;
	m_nasc  = 0;
}

ACatUCAC4::ACatUCAC4(const char *pathdir)
	: ACatalog(pathdir) {
	m_stars = NULL;
	m_asc   = NULL;
	m_nasc  = 0;
}

ACatUCAC4::~ACatUCAC4() {
	if (m_stars) free(m_stars);
	if (m_asc)   free(m_asc);
}

ptr_ucac4_elem ACatUCAC4::GetResult(int &n) {
	n = m_nstars;
	return m_stars;
}

bool ACatUCAC4::LoadAsc() {
	if (m_asc) return true;
	m_nasc  = MILLISEC180 / ucac4_dstep;
	m_nasc *= (MILLISEC360 / ucac4_rstep);
	m_asc = (ptr_ucac4asc) calloc(m_nasc, sizeof(ucac4_asc));
	if (m_asc == NULL) return false;
	char filepath[300];
	sprintf (filepath, "%s/%s/%s", m_pathCat, "u4i", "u4index.asc");
	if (access(filepath, 0)) {
		free (m_asc);
		m_asc = NULL;
		return false;
	}
	// 打开文件
	FILE *fp = fopen(filepath, "r");
	if (fp == NULL) {
		free (m_asc);
		m_asc = NULL;
		return false;
	}
	// 提取文件内容
	char line[40];
	int i = 0;
	unsigned int s, n;
	while (!feof(fp)) {
		if (NULL == fgets(line, 40, fp)) continue;
		sscanf(line, "%u %u", &s, &n);
		m_asc[i].start = s;
		m_asc[i].number= n;
		++i;
	}
	fclose(fp);

	return true;
}

bool ACatUCAC4::FindStar(double ra0, double dec0, double radius) {
	if (!ACatalog::FindStar(ra0, dec0, radius))
		return false;

	LoadAsc();	// 加载快速索引
	vector<ucac4_prime_element> vecrslt;	// 查找到条目的临时缓存区
	int n, i;

	m_csb.zone_seek(ucac4_rstep, ucac4_dstep);

	ra0    *= GtoR;		// 量纲转换, 为后续工作准备
	dec0   *= GtoR;
	radius = radius * GtoR / 60.0;
	// 遍历星表, 查找符合条件的条目
	int zr, zd;		// 赤经赤纬天区编号
	int ZC, ZC0;	// 在索引区中的编号
	double ra, de;	// 星表赤经赤纬
	unsigned int start, number;	// 天区中第一颗星在数据文件中的位置, 和该天区的星数
	FILE *fp;					// 主数据文件访问句柄
	char filepath[300];
	ptr_ucac4_elem buff = NULL;		// 星表数据临时存放地址
	int nelem = 0;					// 星表数据临时存放参考星条目数
	int bytes = (int) sizeof(ucac4_prime_element);

	for (zd = m_csb.zdmin; zd <= m_csb.zdmax; ++zd) {// 遍历赤纬
		sprintf (filepath, "%s/%s/z%03d", m_pathCat, "u4b", zd + 1);
		if ((fp = fopen(filepath, "rb")) == NULL) break;	// 打开文件
		ZC0 = zd * ucac4_zrn;
		for (zr = m_csb.zrmin; zr <= m_csb.zrmax; ++zr) {// 遍历赤经
			ZC = ZC0 + (zr % ucac4_zrn);
			start = m_asc[ZC].start;
			number= m_asc[ZC].number;
			if (number == 0) continue;
			// 为天区数据分配内存
			if (nelem < ((number + 15) & ~15)) {
				free(buff);
				buff = NULL;
				nelem = (number + 15) & ~15;
			}
			if (buff == NULL) buff = (ptr_ucac4_elem) calloc(nelem, sizeof(ucac4_prime_element));
			if (buff == NULL) break;
			// 加载天区数据
			fseek(fp, bytes * start, SEEK_SET);
			fread(buff, bytes, number, fp);
			// 遍历参考星, 检查是否符合查找条件
			for (unsigned int i = 0; i < number; ++i) {
				ra = (double) buff[i].ra / MILLISEC * GtoR;
				de = ((double) buff[i].spd / MILLISEC - 90) * GtoR;
				double v = SphereRange(ra0, dec0, ra, de);
				if (v > radius) continue;
				vecrslt.push_back(*(buff + i));
			}
		}
		fclose(fp);
	}
	if (buff) free(buff);

	// 将查找到的条目存入固定缓存区
	n = (int) vecrslt.size();
	if (AllocBuffer(n)) {
		for (i = 0; i < n; ++i)
			m_stars[i] = vecrslt[i];
	}
	vecrslt.clear();
	return (m_nstars > 0);
}

bool ACatUCAC4::AllocBuffer(int n) {
	if (n > m_max) {// 缓存区不足, 重新分配内存
		m_max = (n + 15) & ~15;	// 条目数为16整数倍
		if (m_stars) {
			free(m_stars);
			m_stars = NULL;
		}
	}
	if ((m_nstars = n) > 0 && m_stars == NULL)
		m_stars = (ptr_ucac4_elem) calloc(m_max, sizeof(ucac4_prime_element));

	return (m_stars != NULL);
}
///////////////////////////////////////////////////////////////////////////////
}
