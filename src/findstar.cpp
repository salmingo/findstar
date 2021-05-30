//============================================================================
// Name        : findstar.cpp
// Author      : Xiaomeng Lu
// Version     :
// Copyright   : Your copyright notice
// Description : 从UCAC4星表中找到符合条件的恒星, 并输出其相关信息
//============================================================================

#include <iostream>
#include <cstdio>
#include <string>
#include <cfitsio/longnam.h>
#include <cfitsio/fitsio.h>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "ADefine.h"
#include "ACatUCAC4.h"
#include "ACatTycho2.h"
#include "ATimeSpace.h"

using namespace std;
using namespace AstroUtil;
using namespace boost::filesystem;
using namespace boost::posix_time;

enum {
	NDX_X,
	NDX_Y,
	NDX_FLUX,
	NDX_MAG,
	NDX_MAGERR,
	NDX_FWHM,
	NDX_ELLIP,
	NDX_BACK,
	NDX_MAX
};

struct param_wcs {
	double x0, y0;	//< XY参考点
	double r0, d0;	//< RA/DEC参考点, 量纲: 弧度
	double cd[2][2];	//< 转换矩阵
	int orderA, orderB;	//< SIP改正阶数
	int ncoefA, ncoefB;	//< SIP系数数量
	double *A, *B;	//< 线性改正系数

public:
	param_wcs() {
		x0 = y0 = 0.0;
		r0 = d0 = 0.0;
		orderA = orderB = 0;
		ncoefA = ncoefB = 0;
		A = B = NULL;
	}

	virtual ~param_wcs() {
		if (A) {
			delete[] A;
			A = NULL;
		}
		if (B) {
			delete[] B;
			B = NULL;
		}
	}

protected:
	/*!
	 * @brief 计算SIP改正模型中与阶数对应的系数数量
	 * @return
	 * 系数数量
	 */
	int term_count(int order) {
		return (order + 1) * (order + 2) / 2;
	}

	/*!
	 * @brief 为SIP改正系数分配存储空间
	 * @param n   系数数量
	 * @param ptr 系数存储地址
	 */
	void alloc_coef(int n, double **ptr) {
		if ((ptr == &A && n != ncoefA) || (ptr == &B && n != ncoefA)) {
			if (ptr == &A)
				ncoefA = n;
			else
				ncoefB = n;
			if (*ptr) {
				delete[] (*ptr);
				(*ptr) = NULL;
			}
		}
		if (*ptr == NULL)
			(*ptr) = new double[n];
	}

	void plane_to_wcs(double xi, double eta, double &ra, double &dec) {
		double fract = cos(d0) - eta * sin(d0);
		ra = cyclemod(r0 + atan2(xi, fract), A2PI);
		dec = atan2(((eta * cos(d0) + sin(d0)) * cos(ra - r0)), fract);
	}

	double poly_val(double x, double y, double *coef, int order) {
		int i, j, k, m;
		double val(0.0), t, px(1.0), py;

		for (i = 0, k = 0; i <= order; ++i) {
			for (j = 0, py = 1.0, t = 0.0, m = order - i; j <= m; ++j, ++k) {
				t += coef[k] * py;
				py *= y;
			}

			val += t * px;
			px *= x;
		}

		return val;
	}

	void project_correct(double &x, double &y) {
		double dx(0.0), dy(0.0);
		dx = poly_val(x, y, A, orderA);
		dy = poly_val(x, y, B, orderB);
		x += dx;
		y += dy;
	}

public:
	bool load_wcs(const string &filepath) {
		fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
		char key[10];
		int status(0), ncoef, i, j, k, m;

		fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRPIX1", &x0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRPIX2", &y0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRVAL1", &r0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CRVAL2", &d0, NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD1_1", &cd[0][0], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD1_2", &cd[0][1], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD2_1", &cd[1][0], NULL, &status);
		fits_read_key(fitsptr, TDOUBLE, "CD2_2", &cd[1][1], NULL, &status);

		fits_read_key(fitsptr, TINT, "A_ORDER", &orderA, NULL, &status);
		if (status)
			return false;
		ncoef = term_count(orderA);
		alloc_coef(ncoef, &A);
		for (i = 0, k = 0; i <= orderA; ++i) {
			for (j = 0, m = orderA - i; j <= m; ++j, ++k) {
				sprintf(key, "A_%d_%d", i, j);
				fits_read_key(fitsptr, TDOUBLE, key, A + k, NULL, &status);
			}
		}

		fits_read_key(fitsptr, TINT, "B_ORDER", &orderB, NULL, &status);
		if (status)
			return false;
		ncoef = term_count(orderB);
		alloc_coef(ncoef, &B);
		for (i = 0, k = 0; i <= orderB; ++i) {
			for (j = 0, m = orderB - i; j <= m; ++j, ++k) {
				sprintf(key, "B_%d_%d", i, j);
				fits_read_key(fitsptr, TDOUBLE, key, B + k, NULL, &status);
			}
		}

		fits_close_file(fitsptr, &status);

		r0 *= D2R;
		d0 *= D2R;
		return !status;
	}

	void image_to_wcs(double x, double y, double &ra, double &dec) {
		double xi, eta;
		x -= x0;
		y -= y0;
		project_correct(x, y);
		xi = (cd[0][0] * x + cd[0][1] * y) * D2R;
		eta = (cd[1][0] * x + cd[1][1] * y) * D2R;
		plane_to_wcs(xi, eta, ra, dec);
		ra *= R2D;
		dec *= R2D;
	}
};

//////////////////////////////////////////////////////////////////////////////
/* 临时以固定数组管理坏像素 */
int bad_col[] = {
	1380
};

int bad_pixel[][2] = {
	{ 943,  179},
	{3568, 1069},
	{ 840, 1201},
	{3976, 1210},
	{2989, 1236},
	{2404, 2307},
	{2458, 2336},
	{1867, 2340},
	{3226, 2894},
	{3227, 2894},
	{3276, 2908},
	{3277, 2908},
	{3319, 2942},
	{3232, 3375},
	{3794, 3408},
	{4051, 3458},
	{4041, 3473},
	{3733, 3800},
	{1509, 3953}
};

/*!
 * @brief 检查是否坏像素
 * @note
 * bad_pixel[][]先按[][1]排序, 若[1]相同, 则按[0]排序
 */
bool is_badpixel(double x, double y) {
	int x0 = int(x + 0.5);
	int y0 = int(y + 0.5);
	int n = sizeof(bad_pixel) / sizeof(int) / 2;
	int low(0), high(n - 1), now;

	if (y0 < bad_pixel[low][1] || y0 > bad_pixel[high][1]) return false;
	while (low < n && bad_pixel[low][1] < y0) ++low;
	while (high < n && bad_pixel[high][1] > y0) --high;
	if (low > high) return false;
	for (now = low; now <= high; ++now) {
		if (x0 == bad_pixel[now][0]) return true;
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @brief 计算大气质量
 * @param z 天顶距, 量纲: 弧度
 * @return
 * 与天顶距对应的大气质量
 */
double airmass(double z) {
	double h = (D2R * 90.0 - z) * R2D; // h: 视高度角, 量纲: 角度
	double am = 1.0 / sin((h + 244.0 / (165.0 + 47.0 * pow(h, 1.1))) * D2R);
	return am;
}

// t: 儒略历元-2000
bool resolve_date_obs(const string &filepath, double &t) {
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int status(0);
	char dateobs[30], timeobs[30];
	string datefull;
	double expt;

	fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
	fits_read_key(fitsptr, TSTRING, "DATE-OBS", dateobs, NULL, &status);
	datefull = dateobs;
	if (datefull.find('T') == string::npos) {
		fits_read_key(fitsptr, TSTRING, "TIME-OBS", timeobs, NULL, &status);
		datefull += "T";
		datefull += timeobs;
	}
	fits_read_key(fitsptr, TDOUBLE, "EXPTIME", &expt, NULL, &status);
	fits_close_file(fitsptr, &status);
	if (status) return false;

	ptime tmobs = from_iso_extended_string(datefull);
	ptime tmmid = tmobs + milliseconds(int(expt * 500));
	t = (tmmid.date().modjulian_day() + tmmid.time_of_day().total_milliseconds() / 86400000.0 - MJD2K) / 365.25;

	return true;
}

/*
int main(int argc, char **argv) {
	path pathroot = argc < 2 ? "./" : argv[2];
	ACatUCAC4 ucac4;
	ucac4item_ptr stars, star, ptr;
	param_wcs wcs;

	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");

	char line[200], ch;
	char seps[] = " ";
	char *token;
	double buff[NDX_MAX];
	double t, ra, dec, racat, decat;
	int i, n, n0(0), n1(0), n2(0);
	FILE *fpstat = fopen("stat.txt", "w");

	for (directory_iterator x = directory_iterator(pathroot); x != directory_iterator(); ++x) {
		if (x->path().extension().string() != ".fit") continue;
		path pathfit = x->path();
		path pathwcs = pathfit;
		path pathcat = pathfit;
		pathwcs.replace_extension(path(".wcs"));
		pathcat.replace_extension(path(".cat"));

		cout << pathfit.filename().string() << "\t\t";
		if (!(exists(pathwcs) && exists(pathcat))) {
			cout << "Absent" << endl;
		}
		else if (!wcs.load_wcs(pathwcs.string())) {
			cout << "WCS Error" << endl;
		}
		else {
			FILE *fpcat = fopen(pathcat.string().c_str(), "r");
			if (!fpcat) {
				cout << "CAT Error" << endl;
			}
			else {
				resolve_date_obs(pathfit.string(), t);
				t *= 1E-4; // 1E-4: 0.1毫角秒=>角秒系数

				path pathmatch = pathfit;
				FILE *fpmatch;

				pathmatch.replace_extension(path(".match"));
				fpmatch = fopen(pathmatch.c_str(), "w");
				n0 = n1 = n2 = 0;
				while (!feof(fpcat)) {
					if (fgets(line, 200, fpcat) == NULL || line[0] == '#') continue;
					n = strlen(line);
					if ((ch = line[n - 1]) == '\r' || ch == '\n') line[n - 1] = 0;
					if ((ch = line[n - 2]) == '\r' || ch == '\n') line[n - 2] = 0;

					i = -1;
					token = strtok(line, seps);
					while (token && ++i < NDX_MAX) {
						buff[i] = atof(token);
						token = strtok(NULL, seps);
					}
					if (abs(int(buff[NDX_X] + 0.5 - 1380)) < 2) continue; // 坏列
					if (is_badpixel(buff[NDX_X], buff[NDX_Y])) continue; // 坏像素
					if (buff[NDX_FLUX] < 1.0 || buff[NDX_FWHM] <= 0.5 || buff[NDX_BACK] > 50000.0) continue;

					++n0;
					wcs.image_to_wcs(buff[NDX_X], buff[NDX_Y], ra, dec);
					if (ucac4.FindStar(ra, dec, 0.5)) {
						stars = ucac4.GetResult(n);
						if (n == 1) {
							++n1;
							star = stars;
						}
						else {
							++n2;
							double dx, dy, dxy2, dxy2min(1E30);
							for (i = 0, ptr = stars; i < n; ++i, ++ptr) {
								decat = (double) star->spd / MILLIAS - 90.0;
								racat = (double) star->ra / MILLIAS + star->pmrac * t * AS2D / cos(decat * D2R);
								decat += (star->pmdc * t * AS2D);

								dx = ra - racat;
								dy = dec - decat;
								if (dx < -180.0) dx += 360.0;
								else if (dx > 180.0) dx -= 360.0;
								dx *= cos(decat * D2R);
								dxy2 = dx * dx + dy * dy;
								if (dxy2 < dxy2min) {
									dxy2min = dxy2;
									star = ptr;
								}
							}
						}
						decat = (double) star->spd / MILLIAS - 90.0;
						racat = (double) star->ra / MILLIAS + star->pmrac * t * AS2D / cos(decat * D2R);
						decat += (star->pmdc * t * AS2D);

						fprintf(fpmatch, "%8.3f %8.3f %10.6f %10.6f %10.6f %10.6f %7.3f %7.3f\r\n",
								buff[NDX_X], buff[NDX_Y], ra, dec, racat, decat,
								buff[NDX_MAG], stars->apasm[3] * 0.001);
					}
				}
				fclose(fpmatch);
				fclose(fpcat);

				printf("%4d  %4d  %4d  %4d\n", n0, n1, n2, n0 - n1 - n2);
				fprintf (fpstat, "%s  %4d  %4d  %4d  %4d\r\n", pathfit.filename().c_str(),
						n0, n1, n2, n0 - n1 - n2);
			}
		}
	}
	fclose(fpstat);

	return 0;
}
*/

///*
int main(int argc, char **argv) {
	if (argc < 3) {
		cout << "************************************************************************" << endl;
		cout << "Usage:" << endl;
		cout << "\t findstar ra dec <radius>" << endl;
		cout << ">> accepted style of R.A.: DDD.DDDD" << endl;
		cout << ">> accepted style of DEC.: ±DD.DDDD" << endl;
		cout << ">> dimension of radius is arc minute, and default value of radius is 1" << endl;
		cout << "************************************************************************" << endl;
		exit(-1);
	}

	double ra0, dec0, radius(1.0);
	double ra, dec;
	int n, i;
	ACatUCAC4 ucac4;
	ACatTycho2 tycho2;
	ATimeSpace ats;

	ra0  = atof(argv[1]);
	dec0 = atof(argv[2]);
	if (argc == 4) radius = atof(argv[3]);
	printf ("UCAC4 result:\n");
	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");
	if (ucac4.FindStar(ra0, dec0, radius)) {
		ucac4item_ptr stars, star;

		printf ("%10s %10s %6s %6s %6s %6s %5s %5s\n",
				"RA_J2000", "DEC_J2000", "pm_ra", "pm_dec",
				"Mag-B", "Mag-V", "Err-B", "Err-V");
		stars = ucac4.GetResult(n);
		for (i = 0, star = stars; i < n; ++i, ++star) {
			ra = (double) star->ra / MILLIAS;
			dec= (double) star->spd / MILLIAS - 90.0;
			printf("%10.6f %10.6f %6.1f %6.1f %6.3f %6.3f %5.2f %5.2f | %10.6f %10.6f\n",
					ra, dec, star->pmrac * 0.1, star->pmdc * 0.1,
					star->apasm[0] * 0.001, star->apasm[1] * 0.001,
					star->apase[0] * 0.01, 	star->apase[1] * 0.01);
		}
	}

	printf ("\nTycho2 result:\n");
	tycho2.SetPathRoot("/Users/lxm/Catalogue/tycho2/tycho2.dat");
	if (tycho2.FindStar(ra0, dec0, radius)) {
		ptr_tycho2_elem stars, star;

		printf ("%10s %10s %6s %6s %6s | %10s %10s\n",
				"RA_J2000", "DEC_J2000", "pm_ra", "pm_dec", "mag");
		stars = tycho2.GetResult(n);
		for (i = 0, star = stars; i < n; ++i, ++star) {
			ra = (double) star->ra / MILLIAS;
			dec= (double) star->spd / MILLIAS - 90.0;
			printf("%10.6f %10.6f %6d %6d %6.3f | %10.6f %10.6f\n",
					ra, dec, star->pmrac, star->pmdc, star->mag * 0.001);
		}
	}

	return 0;
}
// */
