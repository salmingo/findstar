//============================================================================
// Name        : findstar.cpp
// Author      : Xiaomeng Lu
// Version     :
// Copyright   : Your copyright notice
// Description : 从UCAC4星表中找到符合条件的恒星, 并输出其相关信息
//============================================================================

#include <iostream>
#include <vector>
#include <stdio.h>
#include <longnam.h>
#include <fitsio.h>
#include <boost/filesystem.hpp>
#include "ADefine.h"
#include "AMath.h"
#include "ACatUCAC4.h"
#include "ACatTycho2.h"
#include "ATimeSpace.h"

using namespace std;
using namespace AstroUtil;

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

struct file_tag {
	string filename; //< 文件名
	double mjd;	//< 曝光中间时刻对应的修正儒略日
	double lmst;	//< 曝光中间时刻对应的本地平恒星时, 量纲: 弧度
	double expt;	//< 曝光时间
	int wimg, himg;	//< 图像宽度和高度
};

bool inc_mjd(file_tag &x1, file_tag &x2) {
	return (x1.mjd < x2.mjd);
}

typedef vector<file_tag> ftagvec;

bool resolve_date_obs(const string &filepath, ATimeSpace &ats, file_tag &tag) {
	fitsfile *fitsptr;	//< 基于cfitsio接口的文件操作接口
	int status(0);
	char dateobs[40], seps[] = "-T:";
	char *token;
	int year, month, day, hour, minute;
	double second;

	fits_open_file(&fitsptr, filepath.c_str(), 0, &status);
	fits_read_key(fitsptr, TINT, "NAXIS1", &tag.wimg, NULL, &status);
	fits_read_key(fitsptr, TINT, "NAXIS2", &tag.himg, NULL, &status);
	fits_read_key(fitsptr, TSTRING, "DATE-OBS", dateobs, NULL, &status);
	fits_read_key(fitsptr, TDOUBLE, "EXPTIME", &tag.expt, NULL, &status);
	fits_close_file(fitsptr, &status);
	if (status)
		return false;

	token = strtok(dateobs, seps);
	year = atoi(token);
	token = strtok(NULL, seps);
	month = atoi(token);
	token = strtok(NULL, seps);
	day = atoi(token);
	token = strtok(NULL, seps);
	hour = atoi(token);
	token = strtok(NULL, seps);
	minute = atoi(token);
	token = strtok(NULL, seps);
	second = atof(token);

	ats.SetUTC(year, month, day,
			(hour + (minute + (second - tag.expt * 0.5) / 60.0) / 60.0) / 24.0);
	tag.mjd = ats.ModifiedJulianDay();
	tag.lmst = ats.LocalMeanSiderealTime();

	return true;
}

/*!
 * @brief 扫描fits文件, 并按曝光时间增量排序
 * @param pathname 目录名
 * @return
 * 当文件数量过少时返回false
 */
bool scan_fits(const string &pathname, ftagvec &vec, ATimeSpace &ats) {
	namespace fs = boost::filesystem;
	fs::directory_iterator itend = fs::directory_iterator();
	string extdef = ".fit", extname;

	for (fs::directory_iterator x = fs::directory_iterator(pathname);
			x != itend; ++x) {
		extname = x->path().filename().extension().string();
		if (extname.find(extdef) == 0) {
			file_tag tag;
			tag.filename = x->path().filename().string();
			resolve_date_obs(x->path().string(), ats, tag);
			vec.push_back(tag);
		}
	}
	sort(vec.begin(), vec.end(), inc_mjd);
	return vec.size() >= 5;
}

/*
 * 读取目录下的所有wcs文件, 从对应的cat文件中挑选参考星, 并计算
 * 参考星的仪器星等、查找对应的UCAC4星表位置和V星等，计算大气质量.
 * 将所有结果输出至文件output.stars
 */
int main(int argc, char **argv) {
	ATimeSpace ats;
	ats.SetSite(111 + (43.0 + 15.0 / 60.0) / 60.0,
			16 + (27.0 + 4.0 / 60.0) / 60.0, 10.0, 8);

	namespace fs = boost::filesystem;
	ftagvec filevec; // 文件名集合
	fs::path pathroot = argc < 2 ? "./" : argv[2];
	if (!scan_fits(pathroot.string(), filevec, ats)) {
		cout << "require more files to find objects" << endl;
		return -2;
	}

	ACatUCAC4 ucac4;
	ptr_ucac4_elem stars, star;
	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");

	param_wcs wcs;
	char line[200], buff[200];
	char seps[] = " ";
	char *token;
	double x0, y0, dx, dy;
	double x, y, flux, fwhm;
	double ra, dec, rac, decc;
	double azi, ele;
	double mag_inst, mag_vis;
	double range(50.0), min, r;
	int n, i, j, p(0), q(0);
	FILE *output = fopen("output.stars", "w");

	for (ftagvec::iterator it = filevec.begin(); it != filevec.end(); ++it) {
		cout << (*it).filename << endl;
		x0 = (*it).wimg * 0.5 + 0.5;
		y0 = (*it).himg * 0.5 + 0.5;

		fs::path pathfit = pathroot;
		fs::path pathwcs, pathcat;
		pathfit /= (*it).filename;
		pathwcs = pathfit;
		pathcat = pathfit;
		pathwcs.replace_extension(fs::path(".wcs"));
		pathcat.replace_extension(fs::path(".cat"));

		if (!wcs.load_wcs(pathwcs.string())) {
			cout << "failed to load WCS from file: "
					<< pathwcs.filename().string() << endl;
			continue;
		}

		// 4: 将文件名、时间和文件中提取的候选OT赋值给类AFindPV
		FILE *fpcat = fopen(pathcat.string().c_str(), "r");
		if (!fpcat) {
			cout << "failed to open file : " << pathcat.filename().string()
					<< endl;
			continue;
		}

		p = 0;
		while (!feof(fpcat)) {
			if (fgets(line, 200, fpcat) == NULL || line[0] == '#')
				continue;
			n = strlen(line);
			line[n - 1] = 0;
			strcpy(buff, line);

			token = strtok(line, seps);
			x = atof(token);
			token = strtok(NULL, seps);
			y = atof(token);
			token = strtok(NULL, seps);
			flux = atof(token);
			token = strtok(NULL, seps);
			fwhm = atof(token);

//			if (flux < 8000.0 || fwhm < 1.0 || fabs(x - x0) > range
//					|| fabs(y - y0) > range)
//				continue;
			if (flux < 10000.0 || fwhm < 1.0)
				continue;

			wcs.image_to_wcs(x, y, ra, dec);

			if (ucac4.FindStar(ra, dec, 0.2)) {
				stars = ucac4.GetResult(n);
				if (n) {
					min = API;
					ra *= D2R;
					dec *= D2R;
					for (i = 0; i < n; ++i) {
						rac = (double) stars[i].ra / MILLISEC;
						decc = (double) stars[i].spd / MILLISEC - 90.0;
						r = SphereRange(ra, dec, rac * D2R, decc * D2R);
						if (r < min) {
							min = r;
							j = i;
						}
					}
					rac = (double) stars[j].ra / MILLISEC;
					decc = (double) stars[j].spd / MILLISEC - 90.0;

					if (stars[j].apasm[1] < 20000) {
						++p;
						ats.Eq2Horizon((*it).lmst - rac * D2R, decc * D2R, azi,
								ele);
						mag_inst = 25.0 - 2.5 * log10(flux / (*it).expt);
						mag_vis = stars[j].apasm[1] * 0.001;
						// # 星表赤经/赤纬/仪器星等/视星等/星等偏差/天顶距/大气质量
						fprintf(output,
								"%s  %9.5f  %9.5f  %5.2f  %5.2f  %5.2f  %5.2f  %6.3f\r\n",
								buff, rac, decc, mag_inst, mag_vis,
								mag_inst - mag_vis, ele * R2D,
								airmass(D2R * 90.0 - ele));
					}
				}
			}
		}

		fclose(fpcat);
		q += p;
		printf("%d of %d stars being selected\n", p, q);
	}

	fclose(output);
	return 0;
}

/*
 * 读取文件, 从星表中查找对应星的坐标并计算偏差
 */
/*
int main(int argc, char **argv) {
	if (argc < 3) {
		cout << "Usage:" << endl;
		cout << "\t findstar cat_file out_file [exptime = 1.0]" << endl;
		return -1;
	}

	namespace fs = boost::filesystem;
	char line[300], buff[300];
	char *token;
	int pos;
	double x, y, ra0, dec0, ra, dec, r, rmax, flux, fwhm, expt(1.0);
	double azi, ele, mjd, lmst, pmra, pmdec, jy;
	double lgt, lat, alt(10.0);
	double inst_mag, vis_mag;
	int n, i, j;
	ATimeSpace ats;
	ACatUCAC4 ucac4;
	ptr_ucac4_elem stars, star;
	param_wcs wcs;
	FILE *in = fopen(argv[1], "r");
	FILE *ou = fopen(argv[2], "w");
	if (argc == 4)
		expt = atof(argv[3]);

	fs::path pathwcs = argv[1];
	pathwcs.replace_extension(fs::path(".wcs"));
	if (!wcs.load_wcs(pathwcs.string())) {
		cout << "failed to load WCS from " << pathwcs.string() << endl;
		return -2;
	}

	// ...
	lgt = 111 + (43.0 + 15.0 / 60.0) / 60.0;
	lat = 16 + (27.0 + 4.0 / 60.0) / 60.0;
	ats.SetSite(lgt, lat, alt);
	int wimg, himg;
	fs::path pathfit = argv[1];
	pathfit.replace_extension(fs::path(".fit"));
	resolve_date_obs(pathfit.string(), ats, mjd, lmst, expt, wimg, himg);
	jy = (mjd - 51544.5) / 365.25;
	// ...

	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");

	while (!feof(in)) {
		if (fgets(line, 300, in) == NULL)
			continue;
		n = strlen(line);
		line[n - 1] = 0;
		strcpy(buff, line);
		if (line[0] != '#') {
			token = strtok(line, " ");
			x = atof(token);
			token = strtok(NULL, " ");
			y = atof(token);
			token = strtok(NULL, " ");
			flux = atof(token);
			token = strtok(NULL, " ");
			fwhm = atof(token);

			if (flux < 500.0 || fwhm < 1.0)
				continue;

			wcs.image_to_wcs(x, y, ra0, dec0);

			if (ucac4.FindStar(ra0, dec0, 0.1)) {
				stars = ucac4.GetResult(n);
				if (n) {
					rmax = API;
 ra0 *= D2R;
 dec0 *= D2R;
					for (i = 0; i < n; ++i) {
						ra = (double) stars[i].ra / MILLISEC;
						dec = (double) stars[i].spd / MILLISEC - 90.0;
 r = SphereRange(ra0, dec0, ra * D2R, dec * D2R);
						if (r < rmax) {
							rmax = r;
							j = i;
						}
					}

					ra = (double) stars[j].ra / MILLISEC;
					dec = (double) stars[j].spd / MILLISEC - 90.0;
//					pmra = stars[j].pmrac * 1E-4 * jy / 3600.0
 //							/ cos(dec * D2R);
//					pmdec = stars[j].pmdc * 1E-4 * jy;
//					ra += pmra;
//					dec += pmdec;
					if (stars[j].apasm[1] < 20000 && stars[j].apase[1] < 2) {
 ats.Eq2AltAzi(lmst - ra * D2R, dec * D2R, azi, ele);
						inst_mag = 22.0 - 2.5 * log10(flux / expt);
						vis_mag = stars[j].apasm[1] * 0.001;
						// # 星表赤经/赤纬/仪器星等/视星等/星等偏差/天顶距/大气质量
						fprintf(ou,
								"%s  %9.5f  %9.5f  %5.2f  %5.2f  %5.2f  %5.2f  %6.3f\r\n",
								buff, ra, dec, inst_mag, vis_mag,
 inst_mag - vis_mag, 90.0 - ele * R2D,
								airmass(PI90 - ele));
//						fprintf(ou,
//								"%s %9.5f %9.5f %9.5f %9.5f %5.1f %5.1f %5.2f %5.2f\r\n",
 //								buff, ra0 * R2D, dec0 * R2D, ra, dec,
 //								(ra0 * R2D - ra) * 3600.0,
 //								(dec0 * R2D - dec) * 3600.0,
//								25.0 - 2.5 * log10(flux / expt),
//								stars[j].apasm[1] * 0.001);
					}
				}
			}
		}
	}

	fclose(in);
	fclose(ou);

	return 0;
}
 */

/*
 * 命令行输入参数:
 * <1> 赤经. 可接受的格式: HH.HHHH, HH:MM:SS.SS, HHMMSS.SS
 * <2> 赤纬. 可接受的格式: ±DD.DDDD, ±DD:MM:SS.SS, ±DDMMSS.SS
 * <3> 半径, 量纲: 角分. 当缺省时为1
 */
/*
int main(int argc, char **argv) {
	if (argc < 3) {
		cout << "************************************************************************" << endl;
		cout << "Usage:" << endl;
		cout << "\t findstar ra dec <radius>" << endl;
//		cout << ">> accepted style of R.A.: HH.HHHH, HH:MM:SS.SS, HHMMSS.SS" << endl;
//		cout << ">> accepted style of DEC.: ±DD.DDDD, ±DD:MM:SS.S, ±DDMMSS.S" << endl;
		cout << ">> accepted style of R.A.: DDD.DDDD" << endl;
		cout << ">> accepted style of DEC.: ±DD.DDDD" << endl;
		cout << ">> dimension of radius is arc minute, and default value of radius is 1" << endl;
		cout << "************************************************************************" << endl;
		exit(-1);
	}

	double ra0, dec0, radius(1.0);
	double ra, dec;
	ACatUCAC4 ucac4;
	ACatTycho2 tycho2;
	ATimeSpace ats;
	char szRA[20], szDEC[20];

	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");
//	tycho2.SetPathRoot("/Users/lxm/catalogue/tycho2/tycho2.dat");
 //	ats.RAStR2Dbl(argv[1], ra0);
 //	ats.DECStR2Dbl(argv[2], dec0);
	ra0  = atof(argv[1]);
	dec0 = atof(argv[2]);
	if (argc == 4) radius = atof(argv[3]);
	if (ucac4.FindStar(ra0, dec0, radius)) {
//	if (tycho2.FindStar(ra0, dec0, radius)) {
		ptr_ucac4_elem stars, star;
//		ptr_tycho2_elem stars, star;
		int n, i;

		stars = ucac4.GetResult(n);
//		stars = tycho2.GetResult(n);
		for (i = 0, star = stars; i < n; ++i, ++star) {
			ra = (double) star->ra / MILLISEC / 15;
			dec= (double) star->spd / MILLISEC - 90.0;
			ats.RADbl2Str(ra, szRA);
			ats.DECDbl2Str(dec, szDEC);
			printf("%s %s %10.6f %10.6f %6.3f %6.3f %5.2f %5.2f\n",
					szRA, szDEC, ra * 15.0, dec,
					star->apasm[0] * 0.001,
					star->apasm[1] * 0.001, star->apase[0] * 0.01,
					star->apase[1] * 0.01);
		}
	}

	return 0;
}
 */
