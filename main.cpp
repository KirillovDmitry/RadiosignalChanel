/*
Функции main_dll, check, norm_diagram, test предназначены моделирования канала распространения радиосигнала.
Графический интерфейс реализован на LabVEIW, с помощью которого задаются параметры моделирования и передаются
в настоящую dll. Здесь по радиолокационным формулам рассчитываются энергитические характеристики радиосигнала
в каждой точке канала распространения на заданной промежуточной частоте с учетом фазовых соотношений между
переданном и переотраженным сигналом. Dll заполняет переданные массивы значениями амплитуд, фаз принимаемого
сигнала с помощью фазированной антенной решетки.
*/

#include <math.h>
#include <conio.h>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include "energy.h"
#include <complex>
#include "time.h"

using namespace std;

extern "C"
{
_declspec(dllexport) long main_dll(double f, double Fd, double Fdp,
double *A_, double g0_v, double g0_g, double delta, double P,
double *xt, double *yt, double *zt, double *Vxt, double *Vyt, double *Vzt,
double *xr, double *yr, double *zr, double *Vxr, double *Vyr, double *Vzr, double phi0_t,
double wx, double teta0_t, double phi0_r, double teta0_r, double Gt, double Gr, double dteta_t,
double dphi_t, char pz_t, double dteta_r, double dphi_r, char pz_r, double ds, double eps,
double sigma, double n_phi_t, double n_teta_t, double C_phi_t, double C_teta_t, double n_phi_r,
double n_teta_r, double C_phi_r, double C_teta_r, unsigned long N_output, double *dFdop_d,
double *dFdop_r, unsigned long Np, double *A1_, double *P_disp, long it, long ir, double *delay,
double Ct, double Cr, double *F2t, double *F2r, unsigned long n, double *K_d, double *K_r,
double *sin_mean_, double *phi1, double *phi2, long *N_delta_2,
double *A_11, double *A_12, double *A_21, double *A_22, double *A_13, double *A_3, double *R_ves, double *R_AMP,
double *R_FAZA,  double *T_amp, double *T_faza, double d1, double d2,	long size_ves, long size_amp_t,
long size_amp_r, double *PHI_11, double *PHI_12, double *PHI_21, double *PHI_22, double *PHI_13,
double *PHI_3, double *ampl_signal, double *phaza_signal, double EPR, double *T_ves, double *MaxDN_T,
double *MaxDN_R, long number_of_chanel, double *R_amp,double *R_faza, double *Function, 
double *epsilon1, double *epsilon2, double *eps0, double *phi0, double *DEBUG);
_declspec(dllexport) long interpol_faza(double *mas, long size, double* new_mas, long multiply);
_declspec(dllexport) long interpol_amp(double *mas, long size, double* new_mas, long multiply, double *faza);
_declspec(dllexport) long check(double *functional, long size_functional, double phi0_1,
								double phi0_2, double a_1, double a_2, double* phi_1, double* phi_2, double SNR,double *gebug);
_declspec(dllexport) long norm_diagram(double  d1, double  d2, double* T_ves, long T_size_ves, double* R_ves,
									   long R_size_ves, double* T_amp, double* T_faza, long T_size_amp,
									   double* R_amp, double* R_faza, long R_size_amp, double *T_DN_phi, 
									   double *T_DN_eps, double *R_DN_phi, double *R_DN_eps,
									   double *T_Max, double *R_Max, double step, double f, long NumberOfChanel);
_declspec(dllexport) long signal_vxodnoy(long N, double ampl, double Dfrec, double F_d, double T, double t, double *amplituda, double *phasa, double *time);
_declspec(dllexport) long test(double eps1, double eps2, double a1, double a2, double *F); 
}

// основная функция, вычисляющая энергитические характеристика канала распространения радиосигнала
_declspec(dllexport) long main_dll(double f, double Fd, double Fdp,
double *A_, double g0_v, double g0_g, double delta, double P,
double *xt, double *yt, double *zt, double *Vxt, double *Vyt, double *Vzt,
double *xr, double *yr, double *zr, double *Vxr, double *Vyr, double *Vzr, double phi0_t,
double wx, double teta0_t, double phi0_r, double teta0_r, double Gt, double Gr, double dteta_t,
double dphi_t, char pz_t, double dteta_r, double dphi_r, char pz_r, double ds, double eps,
double sigma, double n_phi_t, double n_teta_t, double C_phi_t, double C_teta_t, double n_phi_r,
double n_teta_r, double C_phi_r, double C_teta_r, unsigned long N_output, double *dFdop_d,
double *dFdop_r, unsigned long Np, double *A1_, double *P_disp, long it, long ir, double *delay,
double Ct, double Cr, double *F2t, double *F2r, unsigned long n, double *K_d, double *K_r,
double *sin_mean_, double *PHI1, double *PHI2, long *N_delta_2,
double *A_11, double *A_12, double *A_21, double *A_22, double *A_13, double *A_3, double *R_ves, double *R_AMP,
double *R_FAZA,  double *T_amp, double *T_faza, double d1, double d2,	long size_ves, long size_amp_t, long size_amp_r,
double *PHI_11, double *PHI_12, double *PHI_21, double *PHI_22, double *PHI_13, double *PHI_3, double *ampl_signal,
double *phaza_signal, double EPR, double *T_ves, double *MaxDN_T,double *MaxDN_R, long number_of_chanel, double *R_amp,
double *R_faza, double *Function, double *epsilon1, double *epsilon2, double *EPS0, double *PHI0, double *DEBUG)
{	
/*
описание аргументов функции:
входные параметры:
	f		 	  - частота несущего колебания, Гц
	Fd		 	  - частота дискретизации для прямого и зеркально отражённого сигналов, Гц
	Fdp		 	  - частота дискретизации для рассеянного сигнала, Гц
	Np       	  - отношение частот дискретизации прямого и рассеянного сигналов
	g0_v	 	  - коэффициент обратного рассеяния при нормальном падении для вертикальной поляризации
	g0_g		  - коэффициент обратного рассеяния при нормальном падении для горизонтальной поляризации
	delta	 	  - длина стороны квадрата, на которые разбивается площадь, формирующая рассеянный сигнал, м
	P		 	  - мощность, подводимая к антенне передатчика, Вт
	n		 	  - номер текущей итерации (точки для моделирования)
	xt, yt, zt	  - начальные координаты передатчика в трёхмерной декартовой системе, м
	Vxt, Vyt, Vzt - координаты вектора скорости передатчика в трёхмерной декартовой системе, м/с
	xr, yr, zr	  - начальные координаты приёмника в трёхмерной декартовой системе, м
	Vxr, Vyr, Vzr - координаты вектора скорости приёмника в трёхмерной декартовой системе, м/с
	phi0_t		  - положение оси ДН антенны передатчика, задаваемое углом места в локальной модифицированной сферической системе координат, рад.
	wx			  - угловая скорость вращения передатчика вокруг своей оси (скорость поворота ДН антенны передатчика по углу phi_t), рад/с
	teta0_t		  - положение оси ДН антенны передатчика, задаваемое азимутальным углом в локальной сферической системе координат, рад.
	phi0_r		  - положение оси ДН антенны приёмника, задаваемое азимутальным углом в локальной сферической системе координат, рад.
	teta0_r		  - положение оси ДН антенны приёмника, задаваемое углом места в локальной сферической системе координат, рад.
	Gt			  - коэффициент усиления антенны передатчика (в разах, не в дБ)
	Gr			  - коэффициент усиления антенны приёмника (в разах, не в дБ)
	dteta_t		  - ширина основного лепестка ДН антенны передатчика по азимутальному углу, рад.
	dphi_t		  - ширина основного лепестка ДН антенны передатчика по углу места, рад.
	pz_t		  - поляризация антенны передатчика (относительно оси передатчика)
	dteta_r		  - ширина основного лепестка ДН антенны приёмника в вертикальной плоскости, рад.
	dphi_r		  - ширина основного лепестка ДН антенны приёмника в горизонтальной плоскости, рад.
	pz_r		  - поляризация антенны передатчика (относительно поверхности Земли)
	it			  - изотропность ДН передатчика (1 - изотропная, 0 - нет)
	ir			  - изотропность ДН приёмника (1 - изотропная, 0 - нет)
	Ct			  - параметр, определяющий минимально возможное значение пространственной ДН антенны передатчика
	Cr			  - параметр, определяющий минимально возможное значение пространственной ДН антенны приёмника
	ds			  - среднеквадратичная высота неровностей подстилающей поверхности, м
	eps			  - диэлектрическая проницаемость подстилающей поверхности
	sigma		  - проводимость подстилающей поверхности, См
	n_phi_t		  - степень функции ДН по углу phi (азимутальный угол) антенны передатчика
	n_teta_t	  - степень функции ДН по углу teta (угол места) антенны передатчика
	C_phi_t		  - коэффициент для задания уровня, по которому задаётся ширина ДН по углу phi антенны передатчика
	C_teta_t	  - коэффициент для задания уровня, по которому задаётся ширина ДН по углу teta антенны передатчика
	n_phi_r		  - степень функции ДН по углу phi (азимутальный угол) антенны приёмника
	n_teta_r	  - степень функции ДН по углу teta (угол места) антенны приёмника
	C_phi_r		  - коэффициент для задания уровня, по которому задаётся ширина ДН по углу phi антенны приёмника
	C_teta_r	  - коэффициент для задания уровня, по которому задаётся ширина ДН по углу teta антенны приёмника
	N_output	  - количество точек для моделирования
	
выходные параметры:
	*dFdop_d	  - доплеровский сдвиг прямого сигнала, Гц
	*dFdop_r	  - доплеровский сдвиг квазизеркально отражённого сигнала, Гц
	*P_disp		  - мощность рассеянного сигнала, Вт
	*delay		  - время задержки квазизеркально отражённого сигнала относительно прямого сигнала, сек.
	*F2t		  - значение нормированной ДН антенны передатчика по мощности в направлении прямого луча
	*F2r		  - значение нормированной ДН антенны приёмника по мощности в направлении прямого луча
	*K_d		  - коэффициент рассогласования по амплитуде (из-за несовпадения поляризаций антенны приёмника и приходящего излучения) для прямого сигнала
	*K_r		  - коэффициент рассогласования по амплитуде (из-за несовпадения поляризаций антенны приёмника и приходящего излучения) для квазизеркально отражённого сигнала
	*sin_mean_	  - mean(sin(beta1)+sin(beta2))
	*PHI_11		  - фаза прямого сигнала, рад.
	*PHI_22		  - фаза квазизеркально отражённого сигнала, рад.
	*A_11	 	  - огибающая прямого сигнала РЛС-Ц-РЛС, В (нагрузка 1 Ом)
	*A_22	 	  - огибающая квазизеркально отражённого сигнала РЛС-Ц-РЛС, В (нагрузка 1 Ом)
	*A_12(21)	  - огибающая сигнала,  прямой-квазизеркальный (квазизеркальный-прямой) РЛС-Ц-РЛС, В (нагрузка 1 Ом)
	*A_13    	  - огибающая диффузно-отраженного сигнала РЛС-подстилающая_поверхность-РЛС, В (нагрузка 1 Ом)
	*A_3     	  - огибающая шумового сигнала от подстилающей поверхности на входе антенны РЛС, В (нагрузка 1 Ом)
*/
	unsigned long N = N_output;   // number of samples in output signal
	const double lambda=(3e+8)/f; // длина волны несущего колебания
	long err_no; 				  // код ошибки
	double Td=1/Fd, Tdp=1/Fdp;    // период дискретизации
	
	if(C_phi_t==0 || C_phi_r==0 || C_teta_t==0 || C_teta_r==0) return 4; //ошибка задания параметров ДН
	if(zt[n]==0 || zr[n]==0) return 2; //нулевая высота приёмной или передающей антенны
	int chanel = 0;
	
	struct parameters T,R;    // T - объект структуры, содержащий параметры передатчика, R - объект структуры, содержащий параметры приёмника
	struct com_parameters TR; // структура, содержащая параметры, общие для приемника и передатчика
	
	// заполнение структуры общих параметров
	TR.f=f;
	TR.lambda=3e+8/f;
	TR.Td=1/Fd;
	TR.Tdp=1/Fdp;
	TR.g0_v=g0_v;
	TR.g0_g=g0_g;
	TR.delta=delta;
	TR.eps=eps;
	TR.sigma=sigma;
	TR.ds=ds;
	TR.N=N_output;
	TR.Np=Np;
	TR.EPR = EPR;
	
	// параметры передатчика
	T.G=Gt;
	T.x=xt;
	T.y=yt;
	T.z=zt;
	T.Vx=Vxt;
	T.Vy=Vyt;
	T.Vz=Vzt;
	T.i=it;
	T.pz=pz_t;
	T.teta0=teta0_t;
	T.phi0=phi0_t;
	T.dteta=dteta_t;
	T.dphi=dphi_t;
	T.n_teta=n_teta_t;
	T.n_phi=n_phi_t;
	T.C_teta=C_teta_t;
	T.C_phi=C_phi_t;
	T.C=Ct;

	// параметры приемника
	R.G=Gr;
	R.x=xr;
	R.y=yr;
	R.z=zr;
	R.Vx=Vxr;
	R.Vy=Vyr;
	R.Vz=Vzr;
	R.i=ir;
	R.pz=pz_r;
	R.dteta=dteta_r;
	R.dphi=dphi_r;
	R.n_teta=n_teta_r;
	R.n_phi=n_phi_r;
	R.C_teta=C_teta_r;
	R.C_phi=C_phi_r;
	R.C=Cr;

	dip_vector(&T);//вычисление прямоугольных координат вектора дипольного момента антенны передатчика в системе Oxyz
	dip_vector(&R);//вычисление прямоугольных координат вектора дипольного момента антенны приёмника в системе Oxyz

	// инициализация промежуточных переменных
	double PHI_0 = 0;
	double pfd = 0, Power = 0, Power_r = 0, Power00 = 0, dFdop_d_ = 0, F2t_ = 0, F2r_ = 0, K_d_ = 0, dFdop_d1_ = 0,F2t1_ = 0, F2r1_ = 0, K_d1_ = 0;
	double t_ = 0, t1_ = 0, dt = 0, t2_ =0, dt1 = 0, phi = 0, phi_r = 0, phi_d = 0, phi_d_dT = 0, phi_dT = 0, K_r2_ = 0, phi_r_dT = 0;
	unsigned long n_ = 0, n1_ = 0, n2_ = 0;
	double debug = 0, debug1 = 0, debug2 = 0, D11 = 1, D12 = 1, D21 = 1, D22 = 1, D3 = 1, D13 = 1, D = 1;
	double pfd1 = 0, pfd2 = 0, teta1 = 0, teta2 = 0, phi_ = 0, pfd00 = 0;
	double K_r_, mas_dn0[10] = {0};
	double mas_dn_11[10] = {0}, mas_dn_12[10] = {0}, mas_dn_21[10] = {0}, mas_dn_22[10] = {0}, mas_dn_13[10] = {0};
	double mas_dn_33[10] = {0};

	// моделирование прямого луча (РЛС_Цель):
	err_no=direct_ray_power(&pfd00,&Power00,&phi_d,&phi_d_dT,&t1_,&n1_,&F2t1_,&F2r1_,
		&K_d1_,&dFdop_d1_,&dt1,R,T,TR,n*TR.Td,PHI_0,P, chanel, mas_dn0, &D,  &debug); //прямой луч от РЛС к цели
	
	// моделирование зеркально отражённого луча (РЛС_Цель):
	double A = pfd00*EPR;
	err_no=direct_ray_power(&pfd,&Power,&phi,&phi_dT,&t_,&n_,&F2t_,&F2r_,
		&K_d_,&dFdop_d_,&dt,T,R,TR, t1_ ,phi_d, A, chanel, mas_dn_11, &D11,  &debug);

	if(err_no!=0 && err_no!=10) return err_no;
	if(err_no!=10 && n_<N) //условие того, что задержка прямого луча (выраженная в отсчётах) на данной итерации не превосходит длину массива, в который записываются амплитуды прямого луча
	{
		A_11[n_] = sqrt(2*Power);
		PHI_11[n_] = phi + phi_dT;
	}


	// моделирование луча квазизеркального луча (1-2)
	err_no=reflected_ray_power(&pfd1,&pfd2, &teta1, &teta2,&Power,&phi,&phi_,&phi_dT,
						&t_,&n_,&F2t_,&F2r_,&K_r_,&dFdop_d_,&dt,T,R,TR, t1_, phi_d, A, chanel, mas_dn_12,  &D12, &debug);


	if(err_no!=0 && err_no!=10) return err_no;
	if(err_no!=10 && n_<N) //условие того, что задержка прямого луча (выраженная в отсчётах) на данной итерации не превосходит длину массива, в который записываются амплитуды прямого луча
	{
		A_12[n_] = sqrt(2*Power);
		PHI_12[n_] = phi + phi_dT;
	}
	
	// моделирование луча квазизеркального луча (2-1)
	err_no=reflected_ray_power(&pfd1,&pfd2, &teta1, &teta2,&Power,&phi_r,&phi_,&phi_r_dT,
						&t2_,&n2_,&F2t_,&F2r_,&K_r2_,&dFdop_d_,&dt,R,T,TR,n*TR.Td,PHI_0,P, chanel, mas_dn0, &D,  &debug);

	double B = pfd1*EPR+pfd2*EPR;

	err_no=direct_ray_power(&pfd,&Power,&phi,&phi_dT,&t_,&n_,&F2t_,&F2r_,&K_d_,&dFdop_d_,
		&dt,T,R,TR, t2_, phi_r, B, chanel, mas_dn_21, &D21,  &debug);
		
	if(err_no!=0 && err_no!=10) return err_no;
	if(err_no!=10 && n_<N) 
	{
		A_21[n_] = sqrt(2*Power);
		PHI_21[n_] = phi + phi_dT;
	}

	// моделирование луча квазизеркального луча (2-2)
	err_no=reflected_ray_power(&pfd1,&pfd2, &teta1, &teta2,&Power,&phi,&phi_,&phi_dT,
						&t_,&n_,&F2t_,&F2r_,&K_r_,&dFdop_d_,&dt,T,R,TR, t2_, phi_r, B, chanel, mas_dn_22,  &D22, &debug);

	if(err_no!=0 && err_no!=10) return err_no;
	if(err_no!=10 && n_<N) 
	{
		A_22[n_] = sqrt(2*Power);
		PHI_22[n_] = phi + phi_dT;
	}

	// моделирование рассеянного сигнала (РЛС_Земная поверхность_РЛС)
	err_no=disp_ray_power(0,A_13,sin_mean_, N_delta_2,T,R,TR,t1_,A,n, chanel, mas_dn_13, &D13,  &debug);
	if(n==0) err_no=disp_ray_power(1,A_3,sin_mean_, N_delta_2, R, R, TR, n*TR.Td,P,n, chanel, mas_dn_33, &D3,  &debug);
	
	// моделирование рассеянного сигнала (Земная поверхность_РЛС)
	err_no=disp_ray_power(0,A_3,sin_mean_,N_delta_2,R,T,TR,n*TR.Td,P,n, chanel, mas_dn0, &D, &debug);
	if(err_no!=0 && err_no!=10) return err_no;
	if(err_no!=0 && err_no!=10) return err_no;
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	// методом максимального правдоподобия по принимаемой мощности установливаем направление на источник
	// переотражения зондирующего сигнала
	///////////////////////////////////////////////////////////////////////////////////////////////////

	// вычисленные входные сигналы на антенной решетки преобразуем в комплексный вид
	complex <double> a11(0.0, 0.0), a12(0.0, 0.0), a21(0.0, 0.0), a22(0.0, 0.0), One(0.0, 1.0), One2(0.0, -1.0);
	double Z[10] = {0};
	complex <long double> Y1[5][1] = {One}, Y2[1][5] = {One};

	a11 = (complex <double>)(A_11[n]/sqrt(D11))*exp(One*(complex <double>)PHI_11[n]);
	a12 = (complex <double>)(A_12[n]/sqrt(D12))*exp(One*(complex <double>)PHI_12[n]);
	a21 = (complex <double>)(A_21[n]/sqrt(D21))*exp(One*(complex <double>)PHI_21[n]);
	a22 = (complex <double>)(A_22[n]/sqrt(D22))*exp(One*(complex <double>)PHI_22[n]);
	
	// вычисление смеси комплексного сигнала на входе антенны
	for (int nn = 0; nn<6; nn++)
	{			
		Z[2*nn] = real((a11*(mas_dn_11[2*nn]+One*mas_dn_11[2*nn+1])) + (a12*(mas_dn_12[2*nn]+One*mas_dn_12[2*nn+1]))+ (a21*(mas_dn_21[2*nn]+One*mas_dn_21[2*nn+1]))+ (a22*(mas_dn_22[2*nn]+One*mas_dn_22[2*nn+1])) + (sqrt(A_13[n])/D13*mas_dn_13[2*nn]*Z1_rand)/2);
		Z[2*nn + 1]	= imag((a11*(mas_dn_11[2*nn]+One*mas_dn_11[2*nn+1])) + (a12*(mas_dn_12[2*nn]+One*mas_dn_12[2*nn+1]))+ (a21*(mas_dn_21[2*nn]+One*mas_dn_21[2*nn+1]))+ (a22*(mas_dn_22[2*nn]+One*mas_dn_22[2*nn+1])) + (sqrt(A_13[n])/D13*mas_dn_13[2*nn+1]*Z1_rand)/2);	
	}

	// вычисление угловых координат цели, переотразившей зондирующий сигнал
	double eps1 = 0, eps2 = 0;
	err_no = max_specious_method(&eps1, &eps2, Z, Function, 100, TR);
	if (abs(A_11[n])<=1e-15){
		epsilon1[n] = 0;
		epsilon2[n] = 0;
	}
	else{
		epsilon1[n] =  eps1;
		epsilon2[n] =  eps2;
	}

//***********************************************************
	return err_no; // выход из DLL, возвращаем код ошибки
}


// вспомогательная функция для интерполяции отсчетов амплитуды радиосигнала. работает быстрее встроенной в LabVIEW фунции
_declspec(dllexport) long interpol_amp(double *mas, long size, double* new_mas, long multiply, double *faza)
{
	long i = 0;
	long j = 0;
	double maxx = -1e30;
	long new_size = size*multiply;
	
	for (i = 0; i < size; i++){
		if(abs(mas[i])>maxx)
			maxx = abs(mas[i]);
	}
	for(i = 0; i < size; i++){
		if((faza[i]!= 0)&&((mas[i])!=0)){		
			if ((i<size-1)&&(i>0))
			{

				if (((abs(mas[i])!=0)&&( mas[i+1]!=0)&&( mas[i-1]!=0))&&((abs(faza[i])!=0)&&(faza[i+1]!=0)&&(faza[i-1]!=0))){
					for(j = 0; j<multiply+1; j++)
						new_mas[j + multiply*i] = mas[i] + j*(mas[i+1]- mas[i])/multiply;

				}else if (abs(faza[i-1])==0)
				{
					for(j = 0; j<multiply+1; j++)
						new_mas[j + multiply*i] = mas[i] + j*(mas[i+1]- mas[i])/multiply;

				}else if ( abs(faza[i+1])==0)
				{
					for(j = 1; j<multiply+1; j++)
						new_mas[j + multiply*i] = 0;
				}else if ((abs(mas[i])==0)&&( mas[i+1]==0)&&( mas[i-1]==0))
				{
					for(j = 0; j<multiply+1; j++)
						new_mas[j + multiply*i] = 0;

				}else 
				{
					for(j = 1; j<multiply+1; j++)
						new_mas[j + multiply*i] = 0;
				}

			} else if(i==(size-1)){
				for(j = 0; j<multiply; j++)
					new_mas[i*multiply +j] = mas[i];

			}else if(i==0){
					for(j = 0; j<multiply; j++)
						new_mas[i*multiply +j] = mas[i] + j*(mas[i+1]- mas[i])/multiply;
			}
		} else if((i==0)&& ( abs(faza[i+1])!=0)&&(abs(mas[i+1]-mas[i])<0.5*maxx))
		{
			for(j = 0; j<multiply; j++)
				new_mas[i*multiply +j] = mas[i] + j*(mas[i+1]- mas[i])/multiply;
		} else 
		{
			for(j = 0; j<multiply; j++)
				new_mas[i*multiply +j] = 0;
		}
	}
	return 0;
}

// вспомогательная функция для интерполяции отсчетов фазы радиосигнала. работает быстрее встроенной в LabVIEW фунции
_declspec(dllexport) long interpol_faza(double *mas, long size, double* new_mas, long multiply)
{
	long i = 0;
	long j = 0;
	double maxx = 0;
	long new_size = size*multiply;
	double grad = 0, grad0;
	
	for(i = 0; i < size; i++){
				
			if ((i<size-1)&&(i>0))
			{
				grad0 = (abs(mas[i] - mas[i+1]) + abs(mas[i-1] - mas[i]))/2;
				grad = abs(mas[i] - mas[i+1]);
				if ((abs(mas[i])!=0)&&( mas[i+1]!=0)&&( mas[i-1]!=0)){
					for(j = 0; j<multiply; j++)
						new_mas[j + multiply*i] = mas[i] + j*(mas[i+1]- mas[i])/multiply;

				}else if ((abs(mas[i])==0)&&( mas[i+1]==0)&&( mas[i-1]==0))
				{
					for(j = 0; j<multiply; j++)
						new_mas[j + multiply*i] = 0;

				}else if ((abs(mas[i])!=0)&&( mas[i+1]!=0)&&(grad<=grad0)){
					for(j = 0; j<multiply; j++)
						new_mas[j + multiply*i] = mas[i] + j*(mas[i+1]- mas[i])/multiply;


				}else if (( abs(mas[i])==0)&&( mas[i+1]!=0)&&(grad<grad0)){
					for(j = 0; j<multiply; j++)
						new_mas[j + multiply*i] = mas[i] + j*(mas[i+1]- mas[i])/multiply;

				}else if ( (abs(mas[i])!=0)&&( mas[i+1]==0)&&(grad<grad0)){
					for(j = 0; j<multiply; j++)
						new_mas[j + multiply*i] = mas[i] + j*(mas[i+1]- mas[i])/multiply;
				}else
				{
					for(j = 0; j<multiply; j++)
						new_mas[j + multiply*i] = 0;
				}

			} else if(i==(size-1)){
				for(j = 0; j<multiply; j++)
					new_mas[i*multiply +j] = mas[i];
			}else if(i==0){
				grad = abs(mas[i] - mas[i+1]);
				if (grad<grad0){
					for(j = 0; j<multiply; j++)
					new_mas[i*multiply +j] = mas[i] + j*(mas[i+1]- mas[i])/multiply;
				} else
				{
					for(j = 0; j<multiply; j++)
					new_mas[i*multiply +j] = mas[i];
				
				}
			}
		}

	return 0;
}

// функция вычисляет нормированную диаргамму направленности фазированной антенной решетки
_declspec(dllexport) long norm_diagram(double  d1, double  d2, double* T_ves, long T_size_ves, double* R_ves,
									   long R_size_ves, double* T_amp, double* T_faza, long T_size_amp,
									   double* R_amp, double* R_faza, long R_size_amp, double *T_DN_phi, 
									   double *T_DN_eps, double *R_DN_phi, double *R_DN_eps,
									   double *T_Max, double *R_Max, double step, double f, long NumberOfChanel)
{
/*
 входные параметры:
	f     - несущая частота, Гц
	d1	  - расстояние между вибраторами ФАР в строке, м
	d2	  - расстояние между строками вибраторов в ФАР, м
	T_ves - вектор весовых коэффициетнов вибраторов в строке ФАР в режиме передатчика
	R_ves - вектор весовых коэффициетнов вибраторов в строке ФАР в режиме приемника
	NumberOfChanel - номер угломестного канала
	
 выходные параметры:
	T_DN_phi - диаграмма направленности ФАР в режиме передатчика в азимутальной плоскости
	T_DN_eps - диаграмма направленности ФАР в режиме передатчика в угломестной плоскости
	R_DN_phi - диаграмма направленности ФАР в режиме приемника в азимутальной плоскости
	R_DN_eps - диаграмма направленности ФАР в режиме приемника в угломестной плоскости
	T_Max	 - максимальное значение усиления ФАР в режиме передатчика
	R_Max	 - максимальное значение усиления ФАР в режиме приемника	
*/
	// инициализация переменных для промежуточных вычислений
	step = step*M_PI/180;
	for(int i=0; i < NumberOfChanel; i++){
		R_Max[2*i] = 0;
		R_Max[2*i+1] = 0;
	}
	T_Max[0] = 0;
	T_Max[1] = 0;
	double lamda = (3e+8/f);
	double ves[10] = {0}, amp[10] = {0}, faza[10] = {0};
	for (int i = 0; i < T_size_ves; i++){
		ves[i] = *(T_ves + i);
		amp[i] = *(T_amp + i);
		faza[i] = *(T_faza + i);
	}
	double phi1, phi, phi0, eps1, eps, eps0, Max_phi = 0, Max_eps = 0;
	double left = -M_PI/2, right = M_PI/2;
	int size_ves, count = 0;
	size_ves = T_size_ves;
	int C = 0; 

	complex <double> Sum_phi(0.0, 0.0), Sum_eps(0.0, 0.0);
	complex <double> One(0.0, 1.0);
	complex <double> deg(M_PI/180, 0.0);

	//////////////////////////////////////////////////////////////
	// вычисление ДН в азимутальной плоскости в режиме передатчика
	for(phi = left; phi < right; phi+=step)
		{
			phi1 = phi;
			Sum_phi = (0.0, 0.0);
			for(int i = 0; i < size_ves; i++) 	
			{
				Sum_phi += ves[i]*exp(-2*M_PI*(d1/lamda)*One*((double)i-((double)size_ves-1)/2)*sin(phi1));			
			}
			Sum_phi *= cos((d1/lamda)*phi)*cos((d1/lamda)*phi1);	
					
			T_DN_phi[count] = abs(Sum_phi);
			count++;

			if (Max_phi<abs(Sum_phi)){
				Max_phi = abs(Sum_phi);
				phi0 = phi1;
			}
		
		}
		
	// нахождение максимума ДН в азимутальной плоскости в режиме передатчика	
	C = count;
	left = phi0 - step;
	right = phi0 + step;
	Max_phi = 0;
	for(phi = left; phi<=right; phi+=(step/100))
		{
			phi1 = phi;
			Sum_phi = (0.0, 0.0);
			for(int i = 0; i < size_ves; i++) 	
			{
				Sum_phi += ves[i]*exp(-2*M_PI*(d1/lamda)*One*((double)i-((double)size_ves-1)/2)*sin(phi1));
			}
			Sum_phi *= cos((d1/lamda)*phi1)*cos((d1/lamda)*phi1);	
			if (Max_phi<abs(Sum_phi)){
				Max_phi = abs(Sum_phi);
				phi0 = phi1;
			}
		
		}

	// нормировка ДН к максимальному значению
	*T_Max = Max_phi;
	for( int i = 0; i<count; i++)
		T_DN_phi[i] = T_DN_phi[i]/Max_phi;

	//////////////////////////////////////////////////////////////
	// вычисление ДН в угломестной плоскости в режиме передатчика
	left = -M_PI/2;
	right = M_PI/2;
	count = 0;
	Max_eps = 0;
	
	for(eps = left; eps < right; eps+=step)
		{
			eps1 = eps;
			Sum_eps = (0.0, 0.0);
			for(int n = 0; n < T_size_amp; n++)
			{
				Sum_eps += exp(One*faza[n]*deg)*amp[n]*exp(-2*M_PI*(d2/lamda)*One*((double)n-((double)T_size_amp-1)/2)*sin(eps1));
			}
			Sum_eps *= cos(eps1);
			T_DN_eps[count] = abs(Sum_eps);
			count++;
			if (Max_eps<abs(Sum_eps)){
				Max_eps = abs(Sum_eps);
				eps0 = eps1;
			}	
		}

	// нахождение максимума ДН в угломестной плоскости в режиме передатчика
	left = eps0 - step;
	right = eps0 + step;
	Max_eps = 0;

	for(eps = left; eps<=right; eps+=(step/100))
		{
			eps1 = eps;
			Sum_eps = (0.0, 0.0);
			for(int n = 0; n < T_size_amp; n++)
			{
				Sum_eps += exp(One*faza[n]*deg)*amp[n]*exp(-2*M_PI*(d2/lamda)*One*((double)n-((double)T_size_amp-1)/2)*sin(eps1));
			}
			Sum_eps *= cos(eps1);
			if (Max_eps<abs(Sum_eps)){
				Max_eps = abs(Sum_eps);
				eps0 = eps1;
			}
		}

	// нормировка ДН к максимальному значению
	T_Max[1] = Max_eps;
	for( int i = 0; i<count; i++)
		T_DN_eps[i] = T_DN_eps[i]/Max_eps;

	//////////////////////////////////////////////////////////////
	// для каждого канала вычисляем ДН антенны в режиме приема
	for(int n=0; n < NumberOfChanel; n++){
		for (int i = 0; i < R_size_ves; i++){	
			ves[i] = *(R_ves + i);
		}
	
		for (int i = 0; i < R_size_amp; i++){
			amp[i] = R_amp[(i*NumberOfChanel)+n];
			faza[i] = R_faza[(i*NumberOfChanel)+n];
		}
		left = -M_PI/2;
		right = M_PI/2;
		count = 0;
		Max_eps = 0;
		Max_phi = 0;
		size_ves = R_size_ves;
	
		// вычисление ДН в азимутальной плоскости в режиме приемника
		for(phi = left; phi<=right; phi+=step)
			{
				phi1 = phi ;
				Sum_phi = (0.0, 0.0);
				for(int i = 0; i < size_ves; i++) 	
				{
					Sum_phi += ves[i]*exp(-2*M_PI*(d1/lamda)*One*((double)i-((double)size_ves-1)/2)*sin(phi1));
				}
				Sum_phi *= cos((d1/lamda)*phi1)*cos((d1/lamda)*phi1);	

				R_DN_phi[n*C + count] = abs(Sum_phi);
				count++;

				if (Max_phi<abs(Sum_phi)){
					Max_phi = abs(Sum_phi);
					phi0 = phi1;
				}
			
			}
			
		// нахождение максимума ДН в азимутальной плоскости в режиме приемника
		left = phi0 - step;
		right = phi0 + step;
		Max_phi = 0;
		for(phi = left; phi<=right; phi+=(step/100))
			{
				phi1 = phi ;//+ M_PI/2;
				Sum_phi = (0.0, 0.0);
				for(int i = 0; i < size_ves; i++) 	
				{
					Sum_phi += ves[i]*exp(-2*M_PI*(d1/lamda)*One*((double)i-((double)size_ves-1)/2)*sin(phi1));
				}
				Sum_phi *= cos((d1/lamda)*phi1)*cos((d1/lamda)*phi1);	
				if (Max_phi<abs(Sum_phi)){
					Max_phi = abs(Sum_phi);
					phi0 = phi1;
				}
			
			}
		
		// нормировка ДН
		R_Max[2*n] = Max_phi;
		for( int i = 0; i<count; i++){
			R_DN_phi[n*C+ i] = R_DN_phi[n*C + i]/Max_phi;
		}

		// вычисление ДН в угломестной плоскости в режиме приемника
		left = -M_PI/2;
		right = M_PI/2;
		count = 0;
		Max_eps = 0;
		
		for(eps = left; eps<=right; eps+=step)
			{
				eps1 = eps;
				Sum_eps = (0.0, 0.0);
				for(int nn = 0; nn < R_size_amp; nn++)
				{
					Sum_eps += exp(One*faza[nn]*deg)*amp[nn]*exp(-2*M_PI*(d2/lamda)*One*((double)nn-((double)R_size_amp-1)/2)*sin(eps1));
				}
				Sum_eps *= cos(eps1);

				R_DN_eps[n*C + count] = abs(Sum_eps);
				count++;

				if (Max_eps<abs(Sum_eps)){
					Max_eps = abs(Sum_eps);
					eps0 = eps1;
				}	
			}
			
		// нахождение максимума ДН в угломестной плоскости в режиме приемника
		left = eps0 - step;
		right = eps0 + step;
		Max_eps = 0;

		for(eps = left; eps<=right; eps+=(step/100))
			{
				eps1 = eps ;
				Sum_eps = (0.0, 0.0);
				for(int nn = 0; nn < R_size_amp; nn++)
				{
					Sum_eps += exp(One*faza[nn]*deg)*amp[nn]*exp(-2*M_PI*(d2/lamda)*One*((double)nn-((double)R_size_amp-1)/2)*sin(eps1));
				}
				Sum_eps *= cos(eps1);
				if (Max_eps<abs(Sum_eps)){
					Max_eps = abs(Sum_eps);
					eps0 = eps1;
				}
			
			}
			
		// нормировка ДН
		R_Max[2*n+1] = Max_eps;
			for( int i = 0; i<count; i++)
				R_DN_eps[n*C + i] = R_DN_eps[n*C + i]/Max_eps;
	}

	return 0;
}

// функция генерирования входного сигнала
_declspec(dllexport) long signal_vxodnoy(long N, double ampl, double Dfrec, double F_d, double T, double t, double *amplituda, double *phasa, double *time)
{
/*
 входные параметры:
	N     - количество отсчетов сигнала
	ampl  - амплитуда сигнала
	Dfrec - полоса сигнала
	F_d	  - частота дискретизации сигнала
	t     - текущее время
 выходные параметры:
	amplituda - амплитуда выходного сигнала
	phasa	  - фаза выходного сигнала
	time  	  - массив временных отсчетов

*/
	double Cur_n, NumOfPeriod;
	double Period = T*F_d;
	double Period_Impuls = t*F_d;
	
	long n;
	complex <double> One(0.0, 1.0);

	// для каждого временного отсчета вычисляется амплитуда и фаза выходного сигнала
	for(n = 0; n < N; n++)
	{
		NumOfPeriod = floor(n/(Period));
		Cur_n = n - NumOfPeriod*Period;

		if(Cur_n<=Period_Impuls)
		{
			*(phasa + n) = (M_PI*Dfrec*Cur_n*Cur_n)/(Period_Impuls*F_d) - M_PI*Dfrec*Cur_n/F_d;
			*(amplituda + n) = ampl;
			time[n] = 2*M_PI*Dfrec*Cur_n/Period_Impuls - M_PI*Dfrec;
		}
		else
		{
			*(phasa + n) = 0;
			*(amplituda + n) = 0;
		}

	}
	return 0;
}

// вспомогательная функция, тестирующая метод оценивания угловой координаты лоцируемой цели
_declspec(dllexport)long check(double *functional, long size_functional, double phi0_1, double phi0_2, double a_1, double a_2, double* phi_1, double* phi_2, double SNR,double *gebug)
{
	double Y[10];
	complex <double> Y1[5], Y2[5];
	struct:: com_parameters TR;
	TR.lambda = 0.2;
	TR.d2 = 0.57*(TR.lambda);
	double eps1 = phi0_1/180*M_PI;
	double eps2 = phi0_2/180*M_PI;
	srand((unsigned int)time(NULL));
	complex <double> Sum_phi = (0.0, 0.0), Sum2_phi = (0.0, 0.0), Sum_eps = (0.0, 0.0), Sum2_eps = (0.0, 0.0);
	complex <double> One(0.0, 1.0), One2(0.0, -1.0);
	complex <double> deg(M_PI/180, 0.0);
	double lamda = TR.lambda;
	double P;
	P = 0;
	for(int n = 0; n < 5; n++)
	{	
		Y1[n] = (exp((complex <double>)(-2)*(complex <double>)M_PI*One*(complex <double>)((TR.d2)/(TR.lambda))*(complex <double>)sin(eps1)*(complex <double>)n));//*Sum_phi/MaxDN[2*chanel]);
		Y1[2*n + 1] = imag(exp((complex <double>)(-2)*(complex <double>)M_PI*One*(complex <double>)((TR.d2)/(TR.lambda))*(complex <double>)sin(eps1)*(complex <double>)n));//*Sum_phi/MaxDN[2*chanel]);
	}

	for(int n = 0; n < 5; n++)
	{	
		Y2[n] = (exp((complex <double>)(-2)*(complex <double>)M_PI*One*(complex <double>)((TR.d2)/(TR.lambda))*(complex <double>)sin(eps2)*(complex <double>)n));//*Sum_phi/MaxDN[2*chanel]);
		Y2[2*n + 1] = imag(exp((complex <double>)(-2)*(complex <double>)M_PI*One*(complex <double>)((TR.d2)/(TR.lambda))*(complex <double>)sin(eps2)*(complex <double>)n));//*Sum_phi/MaxDN[2*chanel]);
	}

	for(int n = 0; n < 5; n++)
	{
		Y[2*n] = real(a_1*Y1[n] + a_2*Y2[n] + P*(double)rand()/RAND_MAX);
		Y[2*n + 1] =imag(a_1*Y1[n] + a_2*Y2[n]+ P*(double)rand()/RAND_MAX);
	}

	max_specious_method(phi_1, phi_2, Y, functional, size_functional, TR);

return 0;
}
