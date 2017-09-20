/* $Id: ff.h,v 1.3 2006-09-05 20:46:33 manav Exp $ */

#ifndef	FF_H
#define	FF_H


extern	float	IDilogf( float z[2] );
extern	double	IDilog( double z[2] );
extern	float	Gf( float a, float b, float c, float t );
extern	double	G( double a, double b, double c, double t );
extern	float	Hf( float a, float b, float c, float t );
extern	double	H( double a, double b, double c, double t );
extern	void	Pairf( float ( *c )[2],
		      float p1[3], float p2[3],
		      float q1[3], float q2[3] );
extern	void	Pair( double ( *c )[2],
		      double p1[3], double p2[3],
		      double q1[3], double q2[3] );
extern	float	Areaf( float ( *p )[3], int np );
extern	double	Area( double ( *p )[3], int np );
extern	int	Bilinearf( float ( *c )[2] );
extern	int	Bilinear( double ( *c )[2] );
extern	float	IntegralPlanarf( float ( *c )[2] );
extern	double	IntegralPlanar( double ( *c )[2] );
extern	int	Lcisf( float xax[2], float yax[2],
		      float rad[2], float cnt[2] );
extern	int	Lcis( double xax[2], double yax[2],
		      double rad[2], double cnt[2] );
extern	int	LogSelectf( float rad[2], float cnt[2], float psi[2] );
extern	int	LogSelect( double rad[2], double cnt[2], double psi[2] );
extern	float	RMf( float z[2] );
extern	double	RM( double z[2] );
extern	float	IDilogPathf( int k, float rad[2], float cnt[2],
			    float ( *c )[2] );
extern	double	IDilogPath( int k, double rad[2], double cnt[2],
			    double ( *c )[2] );
extern	float	ILogPartf( int p1, float z[2], float ( *c )[2] );
extern	double	ILogPart( int p1, double z[2], double ( *c )[2] );
extern	float	ILogIntegralf( float ( *c )[2], float s );
extern	double	ILogIntegral( double ( *c )[2], double s );
extern	float	Integralf( float ( *c )[2] );
extern	double	Integral( double ( *c )[2] );
extern	void	FormFactorf( float ( *p )[3], int np,
			    float ( *q )[3], int nq,
			     float* factor12, float* factor21);
extern	void	FormFactor( double ( *p )[3], int np,
			    double ( *q )[3], int nq,
			    double* factor12, double* factor21);

extern	int	fferror;
#define	NO_FF_ERROR	0
#define	INCONSISTENT_K	1
#define	OUT_OF_RANGE_K	2

/*
 * this is really defined in ffp.c but moved
 * out here to support custom calling sequences
 */
#define	ALLCOEFF	21	/* length of c array */

#endif	/* FF_H */
