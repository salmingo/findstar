//============================================================================
// Name        : findstar.cpp
// Author      : Xiaomeng Lu
// Version     :
// Copyright   : Your copyright notice
// Description : 从UCAC4星表中找到符合条件的恒星, 并输出其相关信息
//============================================================================

#include <iostream>
#include <stdio.h>
#include "ACatUCAC4.h"
#include "ACatTycho2.h"
#include "ATimeSpace.h"

using namespace std;
using namespace AstroUtil;

/*
 * 读取文件, 从星表中查找对应星的亮度
 */
int main(int argc, char **argv) {
	if (argc < 3) {
		cout << "Usage:" << endl;
		cout << "\t findstar in_file out_file" << endl;
		return -1;
	}

	char line[300];
	char *token;
	int pos;
	double ra0, dec0;
	int n;
	ACatUCAC4 ucac4;
	ptr_ucac4_elem stars, star;
	FILE *in = fopen(argv[1], "r");
	FILE *ou = fopen(argv[2], "w");

	ucac4.SetPathRoot("/Users/lxm/Catalogue/UCAC4");

	while (!feof(in)) {
		if (fgets(line, 300, in) == NULL) continue;
		n = strlen(line);
		line[n - 1] = 0;
		fprintf(ou, "%s", line);
		if (n < 100) fprintf(ou, "\n");
		else {
			pos = 0;
			token = strtok(line, " ");
			while (token && pos < 14) {
				if (++pos == 13) ra0 = atof(token);
				else if (pos == 14) dec0 = atof(token);
				token = strtok(NULL, " ");
			}

			if (ucac4.FindStar(ra0, dec0, 0.01)) {
				stars = ucac4.GetResult(n);
				if (n == 1) {
					fprintf(ou, "  %.1f  %.1f\n", stars[0].apasm[1] * 0.001, stars[0].apasm[3] * 0.001);
				}
			}
		}
	}

	fclose(in);
	fclose(ou);

	return 0;
}

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
//	ats.RAStr2Dbl(argv[1], ra0);
//	ats.DECStr2Dbl(argv[2], dec0);
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
			printf ("%s %s %10.6f %10.6f %6.3f %6.3f\n",
					szRA, szDEC, ra * 15.0, dec,
					star->apasm[0] * 0.001,
					star->apasm[1] * 0.001);
		}
//			printf ("%s %s %6.3f %6.3f %6.3f %6.3f %6.3f\n",
//					szRA, szDEC,
//					star->apasm[0] * 0.001,
//					star->apasm[1] * 0.001,
//					star->apasm[2] * 0.001,
//					star->apasm[3] * 0.001,
//					star->apasm[4] * 0.001);
//		}
	}

	return 0;
}
*/
