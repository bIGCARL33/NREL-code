
///////////////////////////////////////////////
//         HEADER FILE for SAMPA.C           //
//                                           //
// Solar and Moon Position Algorithm (SAMPA) //
//                   for                     //
//        Solar Radiation Application        //
//                                           //
//              August 1, 2012               //
//                                           //
//   Filename: SAMPA.H                       //
//                                           //
//   Afshin Michael Andreas                  //
//   Afshin.Andreas@NREL.gov (303)384-6383   //
//                                           //
//   Solar Resource and Forecasting Group    //
//   Solar Radiation Research Laboratory     //
//   National Renewable Energy Laboratory    //
//   15013 Denver W Pkwy, Golden, CO 80401   //
//                                           //
//  This code is based on the NREL           //
//  technical report "Solar Eclipse          //
//  Monitoring for Solar Energy Applications //
//  using the Solar and Moon Position 		 //
//  Algorithms" by Ibrahim Reda              //
///////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Usage:                                                               //
//                                                                      //
//   1) Obtain the Solar Position Algorithm from NREL (SPA.C and SPA.H) //
//           http://www.nrel.gov/midc/spa/                              //
//                                                                      //
//   2) In calling program, include this header file,                   //
//      by adding this line to the top of file:                         //
//           #include "sampa.h"                                         //
//                                                                      //
//   3) In calling program, declare the SAMPA structure:                //
//           sampa_data sampa;                                          //
//                                                                      //
//   4) Enter the required input values into SAMPA.SPA structure        //
//      (see below, most input values listed in SPA.H comments)         //
//                                                                      //
//   5) Call the SAMPA calculate function and pass the SAMPA structure  //
//      (prototype is declared at the end of this header file):         //
//           sampa_calculate(&sampa);                                   //
//                                                                      //
//   Output values (listed in comments below) will be computed          //
//   and returned in the passed SAMPA structure.  Output will           //
//   be based on function code selected from enumeration below.         // 
//                                                                      //
//   Note: A non-zero return code from sampa_calculate() indicates that //
//         one of the input values did not pass simple bounds tests.    //
//         The valid input ranges and return error codes are also       //
//         listed in SPA.H comments.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef __sun_and_moon_position_algorithm_header
#define __sun_and_moon_position_algorithm_header

#include "spa.h"

//enumeration for function codes to select desired final outputs from SAMPA
enum 
{
    SAMPA_NO_IRR,  //calculate all values except estimated solar irradiances
    SAMPA_ALL      //calculate all values
};

struct
{
	//-----------------Intermediate MPA OUTPUT VALUES--------------------   
	
	double l_prime;		//moon mean longitude [degrees]
	double d;			//moon mean elongation [degrees]
	double m;			//sun mean anomaly [degrees]
	double m_prime;		//moon mean anomaly [degrees]
	double f;           //moon argument of latitude [degrees]
	double l;			//term l
	double r;			//term r
	double b;			//term b
	double lamda_prime; //moon longitude [degrees]
	double beta;		//moon latitude [degrees]
	double cap_delta;   //distance from earth to moon [kilometers]
	double pi;          //moon equatorial horizontal parallax [degrees] 
	double lamda;       //apparent moon longitude [degrees]
	
    double alpha;       //geocentric moon right ascension [degrees]
    double delta;       //geocentric moon declination [degrees]

    double h;           //observer hour angle [degrees]
    double del_alpha;   //moon right ascension parallax [degrees]
    double delta_prime; //topocentric moon declination [degrees]
    double alpha_prime; //topocentric moon right ascension [degrees]
    double h_prime;     //topocentric local hour angle [degrees]

    double e0;          //topocentric elevation angle (uncorrected) [degrees]
    double del_e;       //atmospheric refraction correction [degrees]
    double e;           //topocentric elevation angle (corrected) [degrees]

    //---------------------Final MPA OUTPUT VALUES------------------------

    double zenith;        //topocentric zenith angle [degrees]
    double azimuth_astro; //topocentric azimuth angle (westward from south) [for astronomers]
    double azimuth;       //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]

} mpa_data; //Moon Position Algorithm (MPA) structure

typedef struct
{
	spa_data spa; //Enter required INPUT VALUES into SPA structure (see SPA.H)
	              //spa.function will be forced to SPA_ZA, therefore slope & azm_rotation not required)
	
	mpa_data mpa; //Moon Position Algorithm structure (defined above)
	
	int function; //Switch to choose functions for desired output (from enumeration)
	
	//---------INPUT VALUES required for estimated solar irradiances--------

	double bird_ozone; //total column ozone thickness [cm] -- range from 0.05 - 0.4
	double bird_pwv;   //total column water vapor [cm] -- range from 0.01 - 6.5 
	double bird_aod;   //broadband aerosol optical depth -- range from 0.02 - 0.5
	double bird_ba;	   //forward scattering factor -- 0.85 recommended for rural aerosols
	double bird_albedo;//ground reflectance -- earth typical is 0.2, snow 0.9, vegitation 0.25 
	
	//---------------------Final SAMPA OUTPUT VALUES------------------------
	
	double ems; //local observed, topocentric, angular distance between sun and moon centers [degrees]
	double rs;	//radius of sun disk [degrees]
	double rm;  //radius of moon disk [degrees]
	
	double a_sul;     //area of sun's unshaded lune (SUL) during eclipse [degrees squared]
	double a_sul_pct; //percent area of SUL during eclipse [percent]
	
	double dni;       //estimated direct normal solar irradiance using SERI/NREL Bird Clear Sky Model [W/m^2]
	double dni_sul;   //estimated direct normal solar irradiance from the sun's unshaded lune [W/m^2]

	double ghi;       //estimated global horizontal solar irradiance using SERI/NREL Bird Clear Sky Model [W/m^2]
	double ghi_sul;   //estimated global horizontal solar irradiance from the sun's unshaded lune [W/m^2]

	double dhi;       //estimated diffuse horizontal solar irradiance using SERI/NREL Bird Clear Sky Model [W/m^2]
	double dhi_sul;   //estimated diffuse horizontal solar irradiance from the sun's unshaded lune [W/m^2]
	
} sampa_data; //Solar and Moon Position Algorithm (SAMPA) structure

//Calculate SAMPA output values (in structure) based on input values passed in structure

#define PI    3.1415926535897932384626433832795028841971

#define COUNT 60

enum {TERM_D, TERM_M, TERM_MPR, TERM_F, TERM_LB, TERM_R, TERM_COUNT};

///////////////////////////////////////////////////////
///  Moon's Periodic Terms for Longitude and Distance
///////////////////////////////////////////////////////
const double ML_TERMS[COUNT][TERM_COUNT]=
{
	{0,0,1,0,6288774,-20905355},
	{2,0,-1,0,1274027,-3699111},
	{2,0,0,0,658314,-2955968},
	{0,0,2,0,213618,-569925},
	{0,1,0,0,-185116,48888},
	{0,0,0,2,-114332,-3149},
	{2,0,-2,0,58793,246158},
	{2,-1,-1,0,57066,-152138},
	{2,0,1,0,53322,-170733},
	{2,-1,0,0,45758,-204586},
	{0,1,-1,0,-40923,-129620},
	{1,0,0,0,-34720,108743},
	{0,1,1,0,-30383,104755},
	{2,0,0,-2,15327,10321},
	{0,0,1,2,-12528,0},
	{0,0,1,-2,10980,79661},
	{4,0,-1,0,10675,-34782},
	{0,0,3,0,10034,-23210},
	{4,0,-2,0,8548,-21636},
	{2,1,-1,0,-7888,24208},
	{2,1,0,0,-6766,30824},
	{1,0,-1,0,-5163,-8379},
	{1,1,0,0,4987,-16675},
	{2,-1,1,0,4036,-12831},
	{2,0,2,0,3994,-10445},
	{4,0,0,0,3861,-11650},
	{2,0,-3,0,3665,14403},
	{0,1,-2,0,-2689,-7003},
	{2,0,-1,2,-2602,0},
	{2,-1,-2,0,2390,10056},
	{1,0,1,0,-2348,6322},
	{2,-2,0,0,2236,-9884},
	{0,1,2,0,-2120,5751},
	{0,2,0,0,-2069,0},
	{2,-2,-1,0,2048,-4950},
	{2,0,1,-2,-1773,4130},
	{2,0,0,2,-1595,0},
	{4,-1,-1,0,1215,-3958},
	{0,0,2,2,-1110,0},
	{3,0,-1,0,-892,3258},
	{2,1,1,0,-810,2616},
	{4,-1,-2,0,759,-1897},
	{0,2,-1,0,-713,-2117},
	{2,2,-1,0,-700,2354},
	{2,1,-2,0,691,0},
	{2,-1,0,-2,596,0},
	{4,0,1,0,549,-1423},
	{0,0,4,0,537,-1117},
	{4,-1,0,0,520,-1571},
	{1,0,-2,0,-487,-1739},
	{2,1,0,-2,-399,0},
	{0,0,2,-2,-381,-4421},
	{1,1,1,0,351,0},
	{3,0,-2,0,-340,0},
	{4,0,-3,0,330,0},
	{2,-1,2,0,327,0},
	{0,2,1,0,-323,1165},
	{1,1,-1,0,299,0},
	{2,0,3,0,294,0},
	{2,0,-1,-2,0,8752}
};
///////////////////////////////////////////////////////
///  Moon's Periodic Terms for Latitude
///////////////////////////////////////////////////////
const double MB_TERMS[COUNT][TERM_COUNT]=
{
	{0,0,0,1,5128122,0},
	{0,0,1,1,280602,0},
	{0,0,1,-1,277693,0},
	{2,0,0,-1,173237,0},
	{2,0,-1,1,55413,0},
	{2,0,-1,-1,46271,0},
	{2,0,0,1,32573,0},
	{0,0,2,1,17198,0},
	{2,0,1,-1,9266,0},
	{0,0,2,-1,8822,0},
	{2,-1,0,-1,8216,0},
	{2,0,-2,-1,4324,0},
	{2,0,1,1,4200,0},
	{2,1,0,-1,-3359,0},
	{2,-1,-1,1,2463,0},
	{2,-1,0,1,2211,0},
	{2,-1,-1,-1,2065,0},
	{0,1,-1,-1,-1870,0},
	{4,0,-1,-1,1828,0},
	{0,1,0,1,-1794,0},
	{0,0,0,3,-1749,0},
	{0,1,-1,1,-1565,0},
	{1,0,0,1,-1491,0},
	{0,1,1,1,-1475,0},
	{0,1,1,-1,-1410,0},
	{0,1,0,-1,-1344,0},
	{1,0,0,-1,-1335,0},
	{0,0,3,1,1107,0},
	{4,0,0,-1,1021,0},
	{4,0,-1,1,833,0},
	{0,0,1,-3,777,0},
	{4,0,-2,1,671,0},
	{2,0,0,-3,607,0},
	{2,0,2,-1,596,0},
	{2,-1,1,-1,491,0},
	{2,0,-2,1,-451,0},
	{0,0,3,-1,439,0},
	{2,0,2,1,422,0},
	{2,0,-3,-1,421,0},
	{2,1,-1,1,-366,0},
	{2,1,0,1,-351,0},
	{4,0,0,1,331,0},
	{2,-1,1,1,315,0},
	{2,-2,0,-1,302,0},
	{0,0,1,3,-283,0},
	{2,1,1,-1,-229,0},
	{1,1,0,-1,223,0},
	{1,1,0,1,223,0},
	{0,1,-2,-1,-220,0},
	{2,1,-1,-1,-220,0},
	{1,0,1,1,-185,0},
	{2,-1,-2,-1,181,0},
	{0,1,2,1,-177,0},
	{4,0,-2,-1,176,0},
	{4,-1,-1,-1,166,0},
	{1,0,1,-1,-164,0},
	{4,0,1,-1,132,0},
	{1,0,-1,-1,-119,0},
	{4,-1,0,-1,115,0},
	{2,-2,0,1,107,0}
};
///////////////////////////////////////////////////////////////////////////////////////////////

double fourth_order_polynomial(double a, double b, double c, double d, double e, double x)
{
    return (((a*x + b)*x + c)*x + d)*x + e;
}

double moon_mean_longitude(double jce)
{
	return limit_degrees(fourth_order_polynomial(
		                 -1.0/65194000, 1.0/538841, -0.0015786, 481267.88123421, 218.3164477, jce));
}

double moon_mean_elongation(double jce)
{
	return limit_degrees(fourth_order_polynomial(
		                 -1.0/113065000, 1.0/545868, -0.0018819, 445267.1114034, 297.8501921, jce));
}

double sun_mean_anomaly(double jce)
{
	return limit_degrees(third_order_polynomial(
		                 1.0/24490000, -0.0001536, 35999.0502909, 357.5291092, jce));
}

double moon_mean_anomaly(double jce)
{
	return limit_degrees(fourth_order_polynomial(
		                 -1.0/14712000, 1.0/69699, 0.0087414, 477198.8675055, 134.9633964, jce));
}

double moon_latitude_argument(double jce)
{
	return limit_degrees(fourth_order_polynomial(
		                 1.0/863310000, -1.0/3526000, -0.0036539, 483202.0175233, 93.2720950, jce));
}

void moon_periodic_term_summation(double d, double m, double m_prime, double f, double jce,
								  const double terms[COUNT][TERM_COUNT], double *sin_sum, double *cos_sum)
{
    int i;
	double e_mult, trig_arg;
	double e  = 1.0 - jce*(0.002516 + jce*0.0000074);

	                  *sin_sum=0;
	if (cos_sum != 0) *cos_sum=0;
    for (i = 0; i < COUNT; i++)
	{
		e_mult   = pow(e, fabs(terms[i][TERM_M]));
		trig_arg = deg2rad(terms[i][TERM_D]*d + terms[i][TERM_M]  *m +
			               terms[i][TERM_F]*f + terms[i][TERM_MPR]*m_prime);
                           *sin_sum += e_mult * terms[i][TERM_LB] *sin(trig_arg);
		if (cos_sum != 0)  *cos_sum += e_mult * terms[i][TERM_R]  *cos(trig_arg);
	}
}

void moon_longitude_and_latitude(double jce, double l_prime, double f, double m_prime, double l, double b,
	                                                                   double *lamda_prime, double *beta)
{
	double a1 = 119.75 +    131.849*jce;
	double a2 =  53.09 + 479264.290*jce;
	double a3 = 313.45 + 481266.484*jce;
	double delta_l =  3958*sin(deg2rad(a1))      + 318*sin(deg2rad(a2))   + 1962*sin(deg2rad(l_prime-f));
	double delta_b = -2235*sin(deg2rad(l_prime)) + 175*sin(deg2rad(a1-f)) +  127*sin(deg2rad(l_prime-m_prime))
		             + 382*sin(deg2rad(a3))      + 175*sin(deg2rad(a1+f)) -  115*sin(deg2rad(l_prime+m_prime));

	*lamda_prime = limit_degrees(l_prime + (l + delta_l)/1000000);
	*beta        = limit_degrees(          (b + delta_b)/1000000);
}

double moon_earth_distance(double r)
{
	return 385000.56 + r/1000;
}

double moon_equatorial_horiz_parallax(double delta)
{
	return rad2deg(asin(6378.14/delta));
}

double apparent_moon_longitude(double lamda_prime, double del_psi)
{
	return lamda_prime + del_psi;
}

double angular_distance_sun_moon(double zen_sun, double azm_sun, double zen_moon, double azm_moon)
{
	double zs = deg2rad(zen_sun);
	double zm = deg2rad(zen_moon);

	return rad2deg(acos(cos(zs)*cos(zm) + sin(zs)*sin(zm)*cos(deg2rad(azm_sun - azm_moon))));
}

double sun_disk_radius(double r)
{
	return 959.63/(3600.0 * r);
}

double moon_disk_radius(double e, double pi, double cap_delta)
{
	return 358473400*(1 + sin(deg2rad(e))*sin(deg2rad(pi)))/(3600.0 * cap_delta);
}

void sul_area(double ems, double rs, double rm, double *a_sul, double *a_sul_pct)
{
	double ems2 = ems*ems;
	double rs2  = rs*rs;
	double rm2  = rm*rm;
	double snum, ai, m, s, h;

	if (ems < (rs + rm))
	{
		if (ems <= fabs(rs - rm))
			ai = PI*rm2;
		else {
			snum =  ems2 + rs2 - rm2;
			m    = (ems2 - rs2 + rm2) / (2*ems);
			s    =              snum  / (2*ems);
			h    = sqrt(4*ems2*rs2 - snum*snum) / (2*ems);
			ai   = (rs2*acos(s/rs) - h * s + rm2*acos(m/rm) - h * m);
		}
	} else ai = 0;

	*a_sul = PI*rs2 - ai;
	if (*a_sul < 0) *a_sul = 0;
	*a_sul_pct = *a_sul * 100.0 / (PI*rs2);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Calculate all MPA parameters and put into structure
// Note: All inputs values (listed in SPA header file) must already be in structure
///////////////////////////////////////////////////////////////////////////////////////////
void mpa_calculate(spa_data *spa, mpa_data *mpa)
{
	mpa->l_prime = moon_mean_longitude(spa->jce);
	mpa->d       = moon_mean_elongation(spa->jce);
	mpa->m       = sun_mean_anomaly(spa->jce);
	mpa->m_prime = moon_mean_anomaly(spa->jce);
	mpa->f       = moon_latitude_argument(spa->jce);

	moon_periodic_term_summation(mpa->d, mpa->m, mpa->m_prime, mpa->f, spa->jce, ML_TERMS, &mpa->l, &mpa->r);
	moon_periodic_term_summation(mpa->d, mpa->m, mpa->m_prime, mpa->f, spa->jce, MB_TERMS, &mpa->b,       0);

	moon_longitude_and_latitude(spa->jce, mpa->l_prime, mpa->f, mpa->m_prime, mpa->l, mpa->b,
		                                                       &mpa->lamda_prime, &mpa->beta);

	mpa->cap_delta = moon_earth_distance(mpa->r);
	mpa->pi = moon_equatorial_horiz_parallax(mpa->cap_delta);

	mpa->lamda    = apparent_moon_longitude(mpa->lamda_prime, spa->del_psi);

    mpa->alpha = geocentric_right_ascension(mpa->lamda, spa->epsilon, mpa->beta);
    mpa->delta = geocentric_declination(mpa->beta, spa->epsilon, mpa->lamda);

	mpa->h  = observer_hour_angle(spa->nu, spa->longitude, mpa->alpha);

    right_ascension_parallax_and_topocentric_dec(spa->latitude, spa->elevation, mpa->pi,
                                mpa->h, mpa->delta, &(mpa->del_alpha), &(mpa->delta_prime));
    mpa->alpha_prime = topocentric_right_ascension(mpa->alpha, mpa->del_alpha);
    mpa->h_prime     = topocentric_local_hour_angle(mpa->h, mpa->del_alpha);

    mpa->e0      = topocentric_elevation_angle(spa->latitude, mpa->delta_prime, mpa->h_prime);
	mpa->del_e   = atmospheric_refraction_correction(spa->pressure, spa->temperature,
                                                         spa->atmos_refract, mpa->e0);
    mpa->e       = topocentric_elevation_angle_corrected(mpa->e0, mpa->del_e);

    mpa->zenith        = topocentric_zenith_angle(mpa->e);
    mpa->azimuth_astro = topocentric_azimuth_angle_astro(mpa->h_prime, spa->latitude, mpa->delta_prime);
    mpa->azimuth       = topocentric_azimuth_angle(mpa->azimuth_astro);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Estimate solar irradiances using the SERI/NREL's Bird Clear Sky Model
///////////////////////////////////////////////////////////////////////////////////////////
void estimate_irr(sampa_data *sampa)
{
	bird_data bird;

	bird.zenith   = sampa->spa.zenith;
	bird.r        = sampa->spa.r;
	bird.pressure = sampa->spa.pressure;
	bird.ozone    = sampa->bird_ozone;
	bird.water    = sampa->bird_pwv;
	bird.taua     = sampa->bird_aod;
	bird.ba       = sampa->bird_ba;
	bird.albedo   = sampa->bird_albedo;
	bird.dni_mod  = sampa->a_sul_pct / 100.0;

	bird_calculate(&bird);

	sampa->dni     = bird.direct_normal;
    sampa->dni_sul = bird.direct_normal_mod;
	sampa->ghi     = bird.global_horiz;
    sampa->ghi_sul = bird.global_horiz_mod;
	sampa->dhi     = bird.diffuse_horiz;
    sampa->dhi_sul = bird.diffuse_horiz_mod;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Calculate all SAMPA parameters and put into structure
// Note: All inputs values (listed in SPA header file) must already be in structure
///////////////////////////////////////////////////////////////////////////////////////////
int sampa_calculate(sampa_data *sampa)
{
    int result;

    sampa->spa.function = SPA_ZA;
    result = spa_calculate(&sampa->spa);

    if (result == 0)
	{
		mpa_calculate(&sampa->spa, &sampa->mpa);

		sampa->ems = angular_distance_sun_moon(sampa->spa.zenith, sampa->spa.azimuth,
			                                   sampa->mpa.zenith, sampa->mpa.azimuth);
		sampa->rs  = sun_disk_radius(sampa->spa.r);
		sampa->rm  = moon_disk_radius(sampa->mpa.e, sampa->mpa.pi, sampa->mpa.cap_delta);

		sul_area(sampa->ems, sampa->rs, sampa->rm, &sampa->a_sul, &sampa->a_sul_pct);

		if (sampa->function == SAMPA_ALL) estimate_irr(sampa);
	}

    return result;
}
///////////////////////////////////////////////////////////////////////////////////////////


#endif