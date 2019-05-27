/*
 * ATimeSpace.cpp
 *
 */
#include "ADefine.h"
#include "ATimeSpace.h"
#include "AMath.h"

namespace AstroUtil
{
///////////////////////////////////////////////////////////////////////////////
ATimeSpace::ATimeSpace()
{
	memset(m_stat, 0, LEN_ATS);
	m_lgt = 0.;
	m_lat = 0.;
	m_alt = 0.;
}

ATimeSpace::~ATimeSpace()
{
}

void ATimeSpace::SetTime(int year, int month, int day, double hour)
{
	double t;

	memset(m_stat, 0, LEN_ATS);
	if (0 == DeltaAT(year, month, day, hour, t)) {
		m_val[DAT] = t;
		m_stat[DAT]= true;
	}
	m_val[JD]  = JulianDay(year, month, day, hour);
	m_val[MJD] = ModifiedJulianDay(m_val[JD]);
	Epoch();
	m_stat[JD] = true;
	m_stat[MJD]= true;
}

void ATimeSpace::SetSite(double lgt, double lat, double alt)
{
	m_lgt = lgt * GtoR;
	m_lat = lat * GtoR;
	m_alt = alt;
	memset(m_stat, 0, LEN_ATS);
}

void ATimeSpace::Hour2HMS(const double hour, int &hh, int &mm, double &ss)
{
	double t = hour;
	hh = (int) t;
	t = (t - hh) * 60;
	mm = (int) t;
	ss = (t - mm) * 60;
}

double ATimeSpace::HMS2Hour(const int hh, const int mm, const double ss)
{
	return (hh + mm / 60.0 + ss/ 3600.0);
}

void ATimeSpace::Deg2DMS(const double deg, int &dd, int &mm, double &ss, int &sign)
{
	double t;
	sign = deg < 0 ? -1 : 1;
	dd = (int) (t = deg * sign);
	t = (t - dd) * 60;
	mm = (int) t;
	ss = (t - mm) * 60;
}

double ATimeSpace::DMS2Deg(const int sign, const int dd, const int mm, const double ss)
{
	return (dd + mm / 60.0 + ss / 3600.0) * sign;
}

bool ATimeSpace::RAStr2Dbl(const char *pszVal, double &val)
{
	if (pszVal == NULL || strlen(pszVal) == 0) return false;
	int len = strlen(pszVal);
	char *buff = new char[len + 1];
	char ch;
	int i, ii;
	double temp;
	bool bDot = false;
	int part = 0;

	for (i = 0, ii = 0, val = 0; i < len; ++i) {// scan input string
		if ((ch = pszVal[i]) >= '0' && ch <= '9') {// is digit
			if (ii == 2 && !bDot) {// scan two sequential digits
				if (part < 2) {
					buff[ii] = 0;
					ii = 0;
					temp = atof(buff);
					if (part == 0 && temp >= 0 && temp < 24) // hour complete
						val = temp;
					else if (part == 1 && temp >= 0) // minute complete
						val += temp / 60.0;
					else return false;

					part++;
				}
				else {
					bDot = true;
					buff[ii++] = '.';
				}
			}
			buff[ii++] = ch;
		}
		else if (ch == ':' || ch ==' ') {// scan separator
			if (bDot || part == 2) // more separator in second section
				return false;
			if (ii != 0) {
				buff[ii] = 0;
				ii = 0;
				temp = atoi(buff);
				if (part == 0 && temp >= 0 && temp < 24) val = temp;
				else if (part == 1 && temp >= 0)         val += temp / 60.0;
				else                                     return false;
			}
			part++;
		}
		else if (ch == '.') {
			if(bDot) break;
			bDot = true;
			buff[ii++] = '.';
		}
		else break;
	}

	buff[ii] = 0;
	if (ii > 1 || (ii == 1 && buff[0] != '.')) {
		temp = atof(buff);
		if (part == 0 && temp >= 0 && temp < 24) val = temp;
		else if (part == 1 && temp >= 0)         val += temp / 60;
		else if (part == 2 && temp >= 0)         val += temp / 3600;
	}
	delete []buff;

	return (part <= 2);
}

bool ATimeSpace::DECStr2Dbl(const char *pszVal, double &val)
{
	if(pszVal == NULL || strlen(pszVal) == 0) return false;
	int len = strlen(pszVal);
	char *buff;
	int i, ii;
	char ch;
	double temp;
	bool bDot = false;
	bool bFlag = false;
	int part = 0;

	i = 0;
	switch (ch = pszVal[0])
	{
	case '+':
		i = 1;
		break;
	case '-':
		bFlag = true;
		i = 1;
		break;
	default:
		if(!isdigit(ch) && ch != '.' && ch != ':' && ch != ' ') return false;
		break;
	}

	buff = new char[len + 1];
	for (ii = 0, temp = 0; i < len; ++i) {
		if ((ch = pszVal[i]) >= '0' && ch <= '9') {
			if (ii == 2 && !bDot) {
				if (part < 2) {
					buff[ii] = 0;
					ii = 0;
					temp = atof(buff);

					if(part == 0)                   val = temp;
					else if (part == 1 && val >= 0) val += temp / 60.0;
					else                            return false;

					part++;
				}
				else {
					bDot = true;
					buff[ii++] = '.';
				}
			}

			buff[ii++] = ch;
		}
		else if (ch == ':' || ch ==' ') {
			if (bDot || part == 2) // more separator in second section
				return false;
			if (ii != 0) {
				buff[ii] = 0;
				ii = 0;
				temp = atoi(buff);
				if (part == 0)                  val = temp;
				else if (part == 1 && val >= 0) val += temp / 60.0;
				else                            return false;
			}
			part++;
		}
		else if (ch == '.') {
			if (bDot) return false;
			bDot = true;
			buff[ii++] = '.';
		}
		else return false;
	}

	buff[ii] = 0;
	if (ii > 1 || (ii == 1 && buff[0] != '.')) {
		temp = atof(buff);
		if(part == 0)                   val = temp;
		else if (part == 1 && val >= 0) val += temp / 60;
		else if (part == 2 && val >= 0) val += temp / 3600;
	}
	if(bFlag) val = -val;

	delete []buff;
	return (part <= 2);
}

void ATimeSpace::RADbl2Str(const double val, char *pszVal)
{
	int hh, mm;
	double ss;
	Hour2HMS(val, hh, mm, ss);
	sprintf(pszVal, "%02d:%02d:%06.3f", hh, mm, ss);
}

void ATimeSpace::DECDbl2Str(const double val, char *pszVal)
{
	int dd, mm, sign;
	double ss;
	Deg2DMS(val, dd, mm, ss, sign);
	if (sign == 1) sprintf(pszVal, "+%02d:%02d:%05.2f", dd, mm, ss);
	else           sprintf(pszVal, "-%02d:%02d:%05.2f", dd, mm, ss);
}

bool ATimeSpace::LeapYear(int year)
{
	return ((year % 400 == 0) || ((year % 4 == 0) && (year % 100 != 0)));
}

int ATimeSpace::DayOfYear(int year, int month, int day)
{
	int K = LeapYear(year) ? 1 : 2;
	int n = (int) (275 * month / 9) - K * (int) ((month + 9) / 12) + day - 30;
	return n;
}

void ATimeSpace::Day2Date(int year, int nth, int &month, int &day)
{
	int K = LeapYear(year) ? 1 : 2;
	if (nth < 32) {
		month = 1;
		day   = nth;
	}
	else {
		month = (int) (9 * (K + nth) / 275 + 0.98);
		day   = nth - (int) (275 * month / 9) + K * (int) ((month + 9) / 12) + 30;
	}
}

int ATimeSpace::DayOfWeek(int year, int month, int day)
{
	int n = (int) (JulianDay(year, month, day, 0) + 1.5);
	return n % 7;
}

double ATimeSpace::Epoch()
{
	if (!m_stat[EPOCH]) {
		m_val[EPOCH] = Epoch(m_val[MJD]);
		m_stat[EPOCH]= true;
	}
	return m_val[EPOCH];
}

double ATimeSpace::Epoch(double mjd)
{
	return 2000 + (mjd - 51544.5) / 365.25;
}

double ATimeSpace::Epoch2JC()
{
	if (!m_stat[JC]) {
		m_val[JC] = (Epoch() - 2000) * 0.01;
		m_stat[JC]= true;
	}
	return m_val[JC];
}

double ATimeSpace::Epoch2JC(double ep)
{
	return (ep - 2000) * 0.01;
}

double ATimeSpace::JulianDay(int year, int month, int day, double hour, int eph)
{
	if (month <= 2) {
		year = year - 1;
		month = month + 12;
	}
	int A = (int) (year / 100);
	int B = 0;
	double jd;

	if (eph == 1)
		B = 2 - A + int(A / 4);
	jd = (int) (365.25 * (year + 4716)) + (int) (30.6001 * (month + 1))
			+ day + hour / 24.0 + B - 1524.5;

	return jd;
}

double ATimeSpace::ModifiedJulianDay(int year, int month, int day, double hour)
{
	return ModifiedJulianDay(JulianDay(year, month, day, hour));
}

double ATimeSpace::ModifiedJulianDay(double jd)
{
	return jd - 2400000.5;
}

void ATimeSpace::JD2YMDH(double jd, int &year, int &month, int &day, double &hour)
{
	int jdn = (int) (jd + 0.5);
	int A, B, C, D, E, t;
	A = jdn;
    if (A >= 2299161) {
    	t = (int) ((jdn - 1867216.25) / 36524.25);
		A = jdn + 1 + t - int(t / 4);
	}
	B = A + 1524;
	C = (int) ((B - 122.1) / 365.25);
	D = (int) (365.25 * C);
	E = (int) ((B - D) / 30.6001);
	day = B - D - (int) (30.6001 * E);
	if (E < 14)
		month = E - 1;
	else
		month = E - 13;
	year = month > 2 ? (C - 4716) : (C - 4715);
	hour = (jd - jdn) * 24;
}

double ATimeSpace::TT2TAI(double tt)
{
	return (tt - TTMTAI / DAYSEC);
}

double ATimeSpace::TAI2TT(double tai)
{
	return (tai + TTMTAI / DAYSEC);
}

int ATimeSpace::DeltaAT(int year, int month, int day, double hour, double &dat)
{
	int iyv = 2015;	// 2015年版本
	/* Reference dates (MJD) and drift rates (s/day), pre leap seconds */
	static const double drift[][2] = {
			{ 37300.0, 0.0012960 },
			{ 37300.0, 0.0012960 },
			{ 37300.0, 0.0012960 },
			{ 37665.0, 0.0011232 },
			{ 37665.0, 0.0011232 },
			{ 38761.0, 0.0012960 },
			{ 38761.0, 0.0012960 },
			{ 38761.0, 0.0012960 },
			{ 38761.0, 0.0012960 },
			{ 38761.0, 0.0012960 },
			{ 38761.0, 0.0012960 },
			{ 38761.0, 0.0012960 },
			{ 39126.0, 0.0025920 },
			{ 39126.0, 0.0025920 }
	};
	int NERA1 = (int) (sizeof(drift) / sizeof(double) / 2);
	/* 日期和时间修正 */
	static const struct {
		int year, month;
		double delat;
	} changes[] = {
			{ 1960,  1,  1.4178180 },
			{ 1961,  1,  1.4228180 },
			{ 1961,  8,  1.3728180 },
			{ 1962,  1,  1.8458580 },
			{ 1963, 11,  1.9458580 },
			{ 1964,  1,  3.2401300 },
			{ 1964,  4,  3.3401300 },
			{ 1964,  9,  3.4401300 },
			{ 1965,  1,  3.5401300 },
			{ 1965,  3,  3.6401300 },
			{ 1965,  7,  3.7401300 },
			{ 1965,  9,  3.8401300 },
			{ 1966,  1,  4.3131700 },
			{ 1968,  2,  4.2131700 },
			{ 1972,  1, 10.0       },
			{ 1972,  7, 11.0       },
			{ 1973,  1, 12.0       },
			{ 1974,  1, 13.0       },
			{ 1975,  1, 14.0       },
			{ 1976,  1, 15.0       },
			{ 1977,  1, 16.0       },
			{ 1978,  1, 17.0       },
			{ 1979,  1, 18.0       },
			{ 1980,  1, 19.0       },
			{ 1981,  7, 20.0       },
			{ 1982,  7, 21.0       },
			{ 1983,  7, 22.0       },
			{ 1985,  7, 23.0       },
			{ 1988,  1, 24.0       },
			{ 1990,  1, 25.0       },
			{ 1991,  1, 26.0       },
			{ 1992,  7, 27.0       },
			{ 1993,  7, 28.0       },
			{ 1994,  7, 29.0       },
			{ 1996,  1, 30.0       },
			{ 1997,  7, 31.0       },
			{ 1999,  1, 32.0       },
			{ 2006,  1, 33.0       },
			{ 2009,  1, 34.0       },
			{ 2012,  7, 35.0       },
			{ 2015,  7, 36.0       }
	};
	const int NDAT = (const int) sizeof(changes) / sizeof(changes[0]);
	int i, j, m;
	double da, mjd;
	dat = da = 0;
	mjd = ModifiedJulianDay(year, month, day, hour);
	if (year < changes[0].year)
		return -1;
	if (year > (iyv + 5))
		j = 1;
	m = 12 * year + month;
	for (i = NDAT - 1; i >= 0; --i)
		if (m >= (12 * changes[i].year + changes[i].month))
			break;
	da = changes[i].delat;
	if (i < NERA1)
		da += (mjd - drift[i][0]) * drift[i][1];
	dat = da;

	return 0;
}

double ATimeSpace::GMST0()
{
	if (!m_stat[GMST]) {
		m_val[GMST] = GMST0(m_val[MJD]);
		m_stat[GMST]= true;
	}
	return m_val[GMST];
}

double ATimeSpace::GMST0(double mjd)
{
	double t = (mjd - 51544.5) / 36525;
	double gmst;

	gmst = reduce(100.46061837 + (36000.770053608 + (0.000387933 - t / 38710000) * t) * t, 360.0);
	return gmst * GtoR;
}

double ATimeSpace::GAST0()
{
	if (!m_stat[GAST]) {
		double gmst = GMST0();
		double nl = NutationLongitude();
		double ee = TrueEclipObliquity();
		double gast = gmst + nl * cos(ee);

		m_val[GAST]  = reduce(gast, PI360);
		m_stat[GAST] = true;
	}

	return m_val[GAST];
}

double ATimeSpace::GAST0(double mjd)
{
	double ep = Epoch(mjd);
	double gmst = GMST0(mjd);
	double nl   = NutationLongitude(ep);
	double ee   = TrueEclipObliquity(ep);
	double gast = gmst + nl * cos(ee);
	return reduce(gast, PI360);
}

double ATimeSpace::LocalMeanSiderialTime()
{
	if (!m_stat[LMST]) {
		m_val[LMST] = LocalMeanSiderialTime(m_val[MJD]);
		m_stat[LMST]= true;
	}
	return m_val[LMST];
}

double ATimeSpace::LocalMeanSiderialTime(double mjd)
{
	double t = (mjd - 51544.5) / 36525;
	double gmst;	// 格林威治平恒星时
	double lmst;	// 本地平恒星时

	gmst = 280.46061837 + (13185000.77005374225 + (0.000387933 - t / 38710000) * t) * t;
	lmst = reduce(gmst * GtoR + m_lgt, PI360);
	return lmst;
}

double ATimeSpace::LocalSiderialTime()
{
	if (!m_stat[LAST]) {
		double nl = NutationLongitude();
		double ee = TrueEclipObliquity();

		m_val[LAST] = reduce(LocalMeanSiderialTime() + nl * cos(ee), PI360);
		m_stat[LAST]= true;
	}
	return m_val[LAST];
}

double ATimeSpace::LocalSiderialTime(double mjd)
{
	double ep = Epoch(mjd);
	double nl = NutationLongitude(ep);
	double ee = TrueEclipObliquity(ep);
	return reduce(LocalMeanSiderialTime(mjd) + nl * cos(ee), PI360);
}

void ATimeSpace::CorrectPM(double ep1, double ra0, double dec0, double pmra, double pmdec, double &ra, double &dec)
{
	double t = Epoch() - ep1;
	double era = t * pmra * 0.001;
	double edec= t * pmdec * 0.001;

	ra = reduce(ra0 + era * StoR, PI360);
	dec= dec0 + edec * StoR;
}

void ATimeSpace::EpochTransform(double ep1, double ra1, double dec1, double &ra2, double &dec2)
{
	double ep2 = Epoch();
	double t1 = (ep1 - 2000) * 0.01;	// 历元1和J2000之间的儒略世纪
	double t2 = (ep2 - ep1)  * 0.01;	// 历元2和历元1之间的儒略世纪

	// 岁差(precession)
	double x = (2306.2181 + 1.39656 * t1 - 0.000139 * pow(t1, 2.0)) * t2
			+ (0.30188 - 0.000344 * t1) * pow(t2, 2.0) + 0.017998 * pow(t2, 3.0);
	double y = (2306.2181 + 1.39656 * t1 - 0.000139 * t1 * t1) * t2
			+ (1.09468 + 0.000066 * t1) * pow(t2, 2.0) + 0.018203 * pow(t2, 3.0);
	double z = (2004.3109 - 0.85330 * t1 - 0.000217 * pow(t1, 2.0)) * t2
			- (0.42665 + 0.000217 * t1) * pow(t2, 2.0) - 0.041833 * pow(t2, 3.0);
	x *= StoR;
	y *= StoR;
	z *= StoR;
	double A = cos(dec1) * sin(ra1 + x);
	double B = cos(z) * cos(dec1) * cos(ra1 + x) - sin(z) * sin(dec1);
	double C = sin(z) * cos(dec1) * cos(ra1 + x) + cos(z) * sin(dec1);
	ra2 = atan2(A, B) + y;
	if (fabs(dec1 - PI90) < GtoR || fabs(dec1 + PI90) < GtoR)
		dec2 = acos(sqrt(A * A + B * B));
	else
		dec2 = asin(C);
	// 章动(nutation)
	ra2 = ra1, dec2 = dec1;
	double lamda, beta;	// 黄道坐标. 使用黄道坐标计算是为了避免赤道坐标在极区的计算错误
	Eq2Eclip(ra2, dec2, lamda, beta);
	lamda += NutationLongitude();
	Eclip2Eq(lamda, beta, ra2, dec2);

	double teo = TrueEclipObliquity();
	double no = NutationObliquity();
	double nl = NutationLongitude();
	printf ("nutation longitude = %f\tnutation latitude = %f\n", nl * RtoS, no * RtoS);
	printf ("eclip obliquity = %f\n", teo * RtoG);
	double era = (cos(teo) + sin(teo) * sin(ra2) * tan(dec2)) * nl - cos(ra2) * tan(dec2) * no;
	double ede = sin(teo) * cos(ra2) * nl + sin(ra2) * no;
	printf ("delta ra = %f\t delta dec = %f\n", era * RtoS, ede * RtoS);
	ra2 = reduce(ra2 + era, PI360);
	dec2= dec2 + ede;
}

double ATimeSpace::MoonLatArg()
{
	if (!m_stat[MLA]) {
		m_val[MLA] = MoonLatArg(m_val[EPOCH]);
		m_stat[MLA]= true;
	}
	return m_val[MLA];
}

double ATimeSpace::MoonLatArg(double ep)
{
	double t = Epoch2JC(ep);
	double v = 93.27191 + 483202.017538 * t - 0.0036825 * pow(t, 2) + pow(t, 3) / 327270;
	return (reduce(v, 360.0) * GtoR);
}

double ATimeSpace::MoonLongAscendingNode()
{
	if (!m_stat[MLAN]) {
		m_val[MLAN] = MoonLongAscendingNode(m_val[EPOCH]);
		m_stat[MLAN]= true;
	}
	return m_val[MLAN];
}

double ATimeSpace::MoonLongAscendingNode(double ep)
{
	double t = Epoch2JC(ep);
	double v = 125.04452 - 1934.136261 * t + 0.0020708 * pow(t, 2) + pow(t, 3) / 450000;
	return (reduce(v, 360) * GtoR);
}

double ATimeSpace::NutationLongitude()
{
	if (!m_stat[NL]) {
		double t = Epoch2JC();
		double D = MoonMeanElongation() * GtoR;
		double M = MeanAnomaly();
		double M1= MoonMeanAnomaly();
		double F = MoonLatArg();
		double ML= MoonLongAscendingNode();

		double v = (-171996 - 174.2 * t) * sin(ML)
		  + ( -13187 -   1.6 * t) * sin(2 * (-D + F + ML))
		  + (  -2274 -   0.2 * t) * sin(2 * (F + ML))
		  + (   2062 +   0.2 * t) * sin(2 * ML)
		  + (   1426 -   3.4 * t) * sin(M)
		  + (    712 +   0.1 * t) * sin(M1)
		  + (   -517 +   1.2 * t) * sin(M + 2 * (-D + F + ML))
		  + (   -386 -   0.4 * t) * sin(2 * F + ML)
		  + (   -301            ) * sin(M1 + 2 * (F + ML))
		  + (    217 -   0.5 * t) * sin(-M + 2 * (-D + F + ML))
		  + (   -158            ) * sin(-2 * D + M1)
		  + (    129 +   0.1 * t) * sin(2 * (-D + F) + ML)
		  + (    123            ) * sin(-M1 + 2 * (F + ML))
		  + (     63            ) * sin(2 * D)
		  + (     63 +   0.1 * t) * sin(M1 + ML)
		  + (    -59            ) * sin(-M1 + 2 * (D + F + ML))
		  + (    -58 -   0.1 * t) * sin(-M1 + ML)
		  + (    -51            ) * sin(M1 + 2 * F + ML)
		  + (     48            ) * sin(2 * (-D + M1))
		  + (     46            ) * sin(2 * (-M1 + F) + ML)
		  + (    -38            ) * sin(2 * (D + F + ML))
		  + (    -31            ) * sin(2 * (M1 + F + ML))
		  + (     29            ) * sin(2 * M1)
		  + (     29            ) * sin(2 * (-D + F + ML) + M1)
		  + (     26            ) * sin(2 * F)
		  + (    -22            ) * sin(2 * (-D + F))
		  + (     21            ) * sin(-M1 + 2 * F + ML)
		  + (     17 -   0.1 * t) * sin(2 * M)
		  + (     16            ) * sin(2 * D - M1 + ML)
		  + (    -16 +   0.1 * t) * sin(2 * (-D + M + F + ML))
		  + (    -15            ) * sin(M + ML)
		  + (    -13            ) * sin(-2 * D + M1 + ML)
		  + (    -12            ) * sin(-M + ML)
		  + (     11            ) * sin(2 * (M1 - F))
		  + (    -10            ) * sin(2 * (D + F) - M1 + ML)
		  + (     -8            ) * sin(2 * (D + F + ML) + M1)
		  + (      7            ) * sin(M + 2 * (F + ML))
		  + (     -7            ) * sin(-2 * D + M + M1)
		  + (     -7            ) * sin(-M + 2 * (F + ML))
		  + (     -7            ) * sin(2 * (D + F) + ML)
		  + (      6            ) * sin(2 * D + M1)
		  + (      6            ) * sin(2 * (-D + M1 + F + ML))
		  + (      6            ) * sin(2 * (-D + F) + M1 + ML)
		  + (     -6            ) * sin(2 * (D - M1) + ML)
		  + (     -6            ) * sin(2 * D + ML)
		  + (      5            ) * sin(-M + M1)
		  + (     -5            ) * sin(2 * (-D + F) - M + ML)
		  + (     -5            ) * sin(-2 * D + ML)
		  + (     -5            ) * sin(2 * (M1 + F) + ML)
		  + (      4            ) * sin(2 * (-D + M1) + ML)
		  + (      4            ) * sin(2 * (-D + F) + M + ML)
		  + (      4            ) * sin(M1 - 2 * F)
		  + (     -4            ) * sin(-D + M1)
		  + (     -4            ) * sin(-2 * D + M)
		  + (     -4            ) * sin(D)
		  + (      3            ) * sin(M1 + 2 * F)
		  + (     -3            ) * sin(2 * (-M1 + F + ML))
		  + (     -3            ) * sin(-D - M + M1)
		  + (     -3            ) * sin(M + M1)
		  + (     -3            ) * sin(-M + M1 + 2 * (F + ML))
		  + (     -3            ) * sin(2 * (D + F + ML) - M - M1)
		  + (     -3            ) * sin(3 * M1 + 2 * (F + ML))
		  + (     -3            ) * sin(2 * (D + F + ML) - M);

		m_val[NL] = v * 0.0001 * GtoR / 3600;
		m_stat[NL]= true;
	}
	return m_val[NL];
}

double ATimeSpace::NutationLongitude(double ep)
{
	double t = Epoch2JC(ep);
	double D = MoonMeanElongation(ep) * GtoR;
	double M = MeanAnomaly(ep);
	double M1= MoonMeanAnomaly(ep);
	double F = MoonLatArg(ep);
	double ML= MoonLongAscendingNode(ep);

	double v = (-171996 - 174.2 * t) * sin(ML)
	  + ( -13187 -   1.6 * t) * sin(2 * (-D + F + ML))
	  + (  -2274 -   0.2 * t) * sin(2 * (F + ML))
	  + (   2062 +   0.2 * t) * sin(2 * ML)
	  + (   1426 -   3.4 * t) * sin(M)
	  + (    712 +   0.1 * t) * sin(M1)
	  + (   -517 +   1.2 * t) * sin(M + 2 * (-D + F + ML))
	  + (   -386 -   0.4 * t) * sin(2 * F + ML)
	  + (   -301            ) * sin(M1 + 2 * (F + ML))
	  + (    217 -   0.5 * t) * sin(-M + 2 * (-D + F + ML))
	  + (   -158            ) * sin(-2 * D + M1)
	  + (    129 +   0.1 * t) * sin(2 * (-D + F) + ML)
	  + (    123            ) * sin(-M1 + 2 * (F + ML))
	  + (     63            ) * sin(2 * D)
	  + (     63 +   0.1 * t) * sin(M1 + ML)
	  + (    -59            ) * sin(-M1 + 2 * (D + F + ML))
	  + (    -58 -   0.1 * t) * sin(-M1 + ML)
	  + (    -51            ) * sin(M1 + 2 * F + ML)
	  + (     48            ) * sin(2 * (-D + M1))
	  + (     46            ) * sin(2 * (-M1 + F) + ML)
	  + (    -38            ) * sin(2 * (D + F + ML))
	  + (    -31            ) * sin(2 * (M1 + F + ML))
	  + (     29            ) * sin(2 * M1)
	  + (     29            ) * sin(2 * (-D + F + ML) + M1)
	  + (     26            ) * sin(2 * F)
	  + (    -22            ) * sin(2 * (-D + F))
	  + (     21            ) * sin(-M1 + 2 * F + ML)
	  + (     17 -   0.1 * t) * sin(2 * M)
	  + (     16            ) * sin(2 * D - M1 + ML)
	  + (    -16 +   0.1 * t) * sin(2 * (-D + M + F + ML))
	  + (    -15            ) * sin(M + ML)
	  + (    -13            ) * sin(-2 * D + M1 + ML)
	  + (    -12            ) * sin(-M + ML)
	  + (     11            ) * sin(2 * (M1 - F))
	  + (    -10            ) * sin(2 * (D + F) - M1 + ML)
	  + (     -8            ) * sin(2 * (D + F + ML) + M1)
	  + (      7            ) * sin(M + 2 * (F + ML))
	  + (     -7            ) * sin(-2 * D + M + M1)
	  + (     -7            ) * sin(-M + 2 * (F + ML))
	  + (     -7            ) * sin(2 * (D + F) + ML)
	  + (      6            ) * sin(2 * D + M1)
	  + (      6            ) * sin(2 * (-D + M1 + F + ML))
	  + (      6            ) * sin(2 * (-D + F) + M1 + ML)
	  + (     -6            ) * sin(2 * (D - M1) + ML)
	  + (     -6            ) * sin(2 * D + ML)
	  + (      5            ) * sin(-M + M1)
	  + (     -5            ) * sin(2 * (-D + F) - M + ML)
	  + (     -5            ) * sin(-2 * D + ML)
	  + (     -5            ) * sin(2 * (M1 + F) + ML)
	  + (      4            ) * sin(2 * (-D + M1) + ML)
	  + (      4            ) * sin(2 * (-D + F) + M + ML)
	  + (      4            ) * sin(M1 - 2 * F)
	  + (     -4            ) * sin(-D + M1)
	  + (     -4            ) * sin(-2 * D + M)
	  + (     -4            ) * sin(D)
	  + (      3            ) * sin(M1 + 2 * F)
	  + (     -3            ) * sin(2 * (-M1 + F + ML))
	  + (     -3            ) * sin(-D - M + M1)
	  + (     -3            ) * sin(M + M1)
	  + (     -3            ) * sin(-M + M1 + 2 * (F + ML))
	  + (     -3            ) * sin(2 * (D + F + ML) - M - M1)
	  + (     -3            ) * sin(3 * M1 + 2 * (F + ML))
	  + (     -3            ) * sin(2 * (D + F + ML) - M);

	return (v * 0.0001 * GtoR / 3600);
}

double ATimeSpace::NutationObliquity()
{
	if (!m_stat[NO]) {
		double t = Epoch2JC();
		double D = MoonMeanElongation() * GtoR;
		double M = MeanAnomaly();
		double M1= MoonMeanAnomaly();
		double F = MoonLatArg();
		double ML= MoonLongAscendingNode();

		double v = (  92025 +   8.9 * t) * cos(ML)
		  + (   5736 -   3.1 * t) * cos(2 * (-D + F + ML))
		  + (    977 -   0.5 * t) * cos(2 * (F + ML))
		  + (   -895 +   0.5 * t) * cos(2 * ML)
		  + (     54 -   0.1 * t) * cos(M)
		  + (     -7            ) * cos(M1)
		  + (    224 -   0.6 * t) * cos(M + 2 * (-D + F + ML))
		  + (    200            ) * cos(2 * F + ML)
		  + (    129 -   0.1 * t) * cos(M1 + 2 * (F + ML))
		  + (    -95 +   0.3 * t) * cos(-M + 2 * (-D + F + ML))
		  + (    -70            ) * cos(2 * (-D + F) + ML)
		  + (    -53            ) * cos(-M1 + 2 * (F + ML))
		  + (    -33            ) * cos(M1 + ML)
		  + (     26            ) * cos(-M1 + 2 * (D + F + ML))
		  + (     32            ) * cos(-M1 + ML)
		  + (     27            ) * cos(M1 + 2 * F + ML)
		  + (    -24            ) * cos(2 * (-M1 + F) + ML)
		  + (     16            ) * cos(2 * (D + F + ML))
		  + (     13            ) * cos(2 * (M1 + F + ML))
		  + (    -12            ) * cos(2 * (-D + F + ML) + M1)
		  + (    -10            ) * cos(-M1 + 2 * F + ML)
		  + (     -8            ) * cos(2 * D - M1 + ML)
		  + (      7            ) * cos(2 * (-D + M + F + ML))
		  + (      9            ) * cos(M + ML)
		  + (      7            ) * cos(-2 * D + M1 + ML)
		  + (      6            ) * cos(-M + ML)
		  + (      5            ) * cos(2 * (D + F) - M1 + ML)
		  + (      3            ) * cos(2 * (D + F + ML) + M1)
		  + (     -3            ) * cos(M + 2 * (F + ML))
		  + (      3            ) * cos(-M + 2 * (F + ML))
		  + (      3            ) * cos(2 * (D + F) + ML)
		  + (     -3            ) * cos(2 * (-D + M1 + F + ML))
		  + (     -3            ) * cos(2 * (-D + F) + M1 + ML)
		  + (      3            ) * cos(2 * (D - M1) + ML)
		  + (      3            ) * cos(2 * D + ML)
		  + (      3            ) * cos(2 * (-D + F) - M + ML)
		  + (      3            ) * cos(-2 * D + ML)
		  + (      3            ) * cos(2 * (M1 + F) + ML);

		m_val[NO] = v * 0.0001 * GtoR / 3600;
		m_stat[NO]= true;
	}
	return m_val[NO];
}

double ATimeSpace::NutationObliquity(double ep)
{
	double t = Epoch2JC(ep);
	double D = MoonMeanElongation(ep) * GtoR;
	double M = MeanAnomaly(ep);
	double M1= MoonMeanAnomaly(ep);
	double F = MoonLatArg(ep);
	double ML= MoonLongAscendingNode(ep);

	double v = (  92025 +   8.9 * t) * cos(ML)
	  + (   5736 -   3.1 * t) * cos(2 * (-D + F + ML))
	  + (    977 -   0.5 * t) * cos(2 * (F + ML))
	  + (   -895 +   0.5 * t) * cos(2 * ML)
	  + (     54 -   0.1 * t) * cos(M)
	  + (     -7            ) * cos(M1)
	  + (    224 -   0.6 * t) * cos(M + 2 * (-D + F + ML))
	  + (    200            ) * cos(2 * F + ML)
	  + (    129 -   0.1 * t) * cos(M1 + 2 * (F + ML))
	  + (    -95 +   0.3 * t) * cos(-M + 2 * (-D + F + ML))
	  + (    -70            ) * cos(2 * (-D + F) + ML)
	  + (    -53            ) * cos(-M1 + 2 * (F + ML))
	  + (    -33            ) * cos(M1 + ML)
	  + (     26            ) * cos(-M1 + 2 * (D + F + ML))
	  + (     32            ) * cos(-M1 + ML)
	  + (     27            ) * cos(M1 + 2 * F + ML)
	  + (    -24            ) * cos(2 * (-M1 + F) + ML)
	  + (     16            ) * cos(2 * (D + F + ML))
	  + (     13            ) * cos(2 * (M1 + F + ML))
	  + (    -12            ) * cos(2 * (-D + F + ML) + M1)
	  + (    -10            ) * cos(-M1 + 2 * F + ML)
	  + (     -8            ) * cos(2 * D - M1 + ML)
	  + (      7            ) * cos(2 * (-D + M + F + ML))
	  + (      9            ) * cos(M + ML)
	  + (      7            ) * cos(-2 * D + M1 + ML)
	  + (      6            ) * cos(-M + ML)
	  + (      5            ) * cos(2 * (D + F) - M1 + ML)
	  + (      3            ) * cos(2 * (D + F + ML) + M1)
	  + (     -3            ) * cos(M + 2 * (F + ML))
	  + (      3            ) * cos(-M + 2 * (F + ML))
	  + (      3            ) * cos(2 * (D + F) + ML)
	  + (     -3            ) * cos(2 * (-D + M1 + F + ML))
	  + (     -3            ) * cos(2 * (-D + F) + M1 + ML)
	  + (      3            ) * cos(2 * (D - M1) + ML)
	  + (      3            ) * cos(2 * D + ML)
	  + (      3            ) * cos(2 * (-D + F) - M + ML)
	  + (      3            ) * cos(-2 * D + ML)
	  + (      3            ) * cos(2 * (M1 + F) + ML);

	return (v * 0.0001 * GtoR / 3600);
}

double ATimeSpace::MeanEclipObliquity()
{
	if (!m_stat[MEO]) {
		m_val[MEO] = MeanEclipObliquity(m_val[EPOCH]);
		m_stat[MEO]= true;
	}
    return m_val[MEO];
}

double ATimeSpace::MeanEclipObliquity(double ep)
{
	double jc = Epoch2JC(ep);	// 儒略世纪
	double u  = jc * 0.01;
	// v的量纲是角秒
	double v = 84381.448 - 4680.93 * u
						 -    1.55 * pow(u, 2.0)
						 + 1999.25 * pow(u, 3.0)
						 -   51.38 * pow(u, 4.0)
						 -  249.67 * pow(u, 5.0)
						 -   39.05 * pow(u, 6.0)
						 +    7.12 * pow(u, 7.0)
						 +   27.87 * pow(u, 8.0)
						 +    5.79 * pow(u, 9.0)
						 +    2.45 * pow(u, 10.0);
	return v * GtoR / 3600;
}

double ATimeSpace::TrueEclipObliquity()
{
	if (!m_stat[TEO]) {
		double meo = MeanEclipObliquity();
		double no  = NutationObliquity();

		m_val[TEO] = reduce(meo + no, PI360);
		m_stat[TEO]= true;
	}
	return m_val[TEO];
}

double ATimeSpace::TrueEclipObliquity(double ep)
{
	double meo = MeanEclipObliquity(ep);
	double no  = NutationObliquity(ep);

	return reduce(meo + no, PI360);
}

double ATimeSpace::MoonMeanElongation()
{
	if (!m_stat[MME]) {
		double t = Epoch2JC();
		double v = 297.85036 + 445267.111480 * t - 0.0019142 * pow(t, 2.0) + pow(t, 3.0) / 189474;
		m_val[MME] = reduce(v, 360.0);
		m_stat[MME]= true;
	}
	return m_val[MME];
}

double ATimeSpace::MoonMeanElongation(double ep)
{
	double t = Epoch2JC(ep);
	double v = 297.85036 + 445267.111480 * t - 0.0019142 * pow(t, 2.0) + pow(t, 3.0) / 189474;
	return reduce(v, 360.0);
}

double ATimeSpace::MeanAnomaly()
{
	if (!m_stat[MA]) {
		m_val[MA] = MeanAnomaly(m_val[EPOCH]);
		m_stat[MA]= true;
	}
	return m_val[MA];
}

double ATimeSpace::MeanAnomaly(double ep)
{
	double t = Epoch2JC(ep);
	double v = 357.52772 + 35999.050340 * t - 0.0001603 * pow(t, 2.0) - pow(t, 3.0) / 300000;
	return (reduce(v, 360) * GtoR);
}

double ATimeSpace::MoonMeanAnomaly()
{
	if (!m_stat[MMA]) {
		m_val[MMA] = MoonMeanAnomaly(m_val[EPOCH]);
		m_stat[MMA]= true;
	}
	return m_val[MMA];
}

double ATimeSpace::MoonMeanAnomaly(double ep)
{
	double t = Epoch2JC(ep);
	double v = 134.96298 + 477198.867398 * t + 0.0086972 * pow(t, 2.0) + pow(t, 3.0) / 56250;
	return reduce(v * GtoR, PI360);
}

void ATimeSpace::Eclip2Eq(double l, double b, double &ra, double &dec)
{
	double eo = TrueEclipObliquity();
	ra = atan2((sin(l) * cos(eo) - tan(b) * sin(eo)), cos(l));
	dec= asin(sin(b) * cos(eo) + cos(b) * sin(eo) * sin(l));
	ra = reduce(ra, PI360);
}

void ATimeSpace::Eq2Eclip(double ra, double dec, double &l, double &b)
{
	double eo = TrueEclipObliquity();
	l = atan2(sin(ra) * cos(eo) + tan(dec) * sin(eo), cos(ra));
	b = asin(sin(dec) * cos(eo) - cos(dec) * sin(eo) * sin(ra));
	l = reduce(l, PI360);
}

void ATimeSpace::AltAzi2Eq(double azi, double alt, double &ha, double &dec)
{
	ha  = atan2(sin(azi), cos(azi) * sin(m_lat) + tan(alt) * cos(m_lat));
	dec = asin(sin(m_lat) * sin(alt) - cos(m_lat) * cos(alt) * cos(azi));
}

void ATimeSpace::Eq2AltAzi(double ha, double dec, double &azi, double &alt)
{
	azi = atan2(sin(ha), cos(ha) * sin(m_lat) - tan(dec) * cos(m_lat));
	alt = asin(sin(m_lat) * sin(dec) + cos(m_lat) * cos(dec) * cos(ha));
	azi = reduce(azi, PI360);
}

double ATimeSpace::ParAngle(double ha, double dec, char mode)
{
	double y = sin(ha);
	double x = tan(m_lat) * cos(dec) - sin(dec) * cos(ha);
	double q = atan2(y, x);
	double azi(0), alt(0);
	if(mode != 'N') {
		Eq2AltAzi(ha, dec, azi, alt);
		if(mode == '+') q += alt;
		else if(mode == '-') q -= alt;
	}

	return reduce(q, PI360);
}

///////////////////////////////////////////////////////////////////////////////
}
