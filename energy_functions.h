#define M_PI 3.14159265358979323846	//определение числа пи
#define ROU(x) (((x)-floor((x))<0.5)?(floor((x))):(ceil((x))))	//определение операции округления к ближайшему целому
#define pow2(x) (x)*(x) //определение операции возведения в квадрат

void dip_vector(struct parameters *S);

//структура параметров, имеющих одинаковые значения для приёмника и передатчика:
struct com_parameters
	{ // определение параметров данной структуры см. в описании аргументов функции main_dll
		double f, lambda, Td, Tdp, g0_v, g0_g, delta, eps, sigma, ds;
		unsigned long N, Np;
		double d1, d2, EPR;
		
	};

//структура параметров, имеющих разные значения для приёмника и передатчика:
struct parameters
	{
		double G; //коэффициент усиления антенны
		double *x, *y, *z; // прямоугольные координаты
		double *Vx, *Vy, *Vz; // координаты вектора скорости
		long i; // изотропность антенны
		char pz; // поляризация антенны
		double xp, yp, zp; // прямоугольные координаты вектора дипольного момента антенны
		double teta0, phi0; // угловые координаты оси диаграммы направленности антенны
		double dteta, dphi; // полуширины диаграммы направленности антенны в двух плоскостях
		double n_teta, n_phi; // степени функций диаграммы направленности антенны по углам teta и phi
		double C_teta, C_phi, C; // коэффициенты для задания уровня, по которому задаётся ширина ДН по углам
								 // teta и phi антенны, и параметр, определяющий минимально возможное значение
								 // пространственной ДН антенны (параметр C)
		double *T_ves, *R_ves, *T_amp, *T_faza, *R_amp, *R_faza, *MaxDN_R, *MaxDN_T;
		long size_ves,  T_size_amp, R_size_amp;
	};


//структура координат приёмника и передатчика:
struct coords
	{
		double xt,yt,zt,xr,yr,zr;
	};

// Функция вычисления прямоугольных координат вектора дипольного момента антенны (дипольное приближение антенны)
void dip_vector(struct parameters *S);

// определение направления вектора E линейно поляризованной волны антенны передатчика в направлении на приёмник
double pn(double xp, double yp, double zp, double r0x, double r0y, double r0z, double *xE, double *yE, double *zE);

// Функция для вычисления доплеровского сдвига частоты в прямом и зеркально отражённом лучах
double doppler(char dr, double t_nTd, struct coords C, struct parameters T, struct parameters R, 
					struct com_parameters TR, double *d_2, double *Ld, double *Lr);
					
// функция вычисляет модуль и аргумент комплексного коэффициента отражения Френеля для TE- или TM-компоненты поля
double abs_z3(double *teta, char pz, double eps, double sigma, double f, double sb, double cb);

// функция вычисляят характеристики сигнала, распространяющиего по направлению прямой видимости
int direct_ray_power(double *pfd,double *P,double *phi,double *phi_dT,double *t_,unsigned long *n_,double *F2t,
					 double *F2r,double *K_d,double *dFdop_d,double *dt,struct parameters T,struct parameters R,
					 struct com_parameters TR,double t_nTd,double phi0,double P0, int chanel, double *mas,
					 double *d, double *debug);
					 
// функция вычисляят характеристики сигнала, распространяющиего через точку квазизеркального отражения
int reflected_ray_power(double *pfd1,double *pfd2,double *teta1,double *teta2,double *P,double *phi,double *phi_,
						double *phi_dT,double *t_,unsigned long *n_,double *F2t,double *F2r,double *K_r,
						double *dFdop_r,double *dt,struct parameters T,struct parameters R,struct com_parameters TR,
						double t_nTd,double phi0,double P0, int chanel, double *mas, double *d, double *debug);

// функция вычисляят характеристики диффузно переотраженного сигнала
int disp_ray_power(int s,double *P_disp,double *sin_mean_,long *N_delta_2,struct parameters T,
				   struct parameters R,struct com_parameters TR,double t_nTd,double P0,unsigned long n, int chanel,
				   double *mas, double *d, double *debug);

// двумерная и одномерная функции для вычисления направления на цель методом максимального правдоподобия
int max_specious_method(double *eps1, double *eps2, double *Y1, double *Y2, long size,struct com_parameters TR);
int max_specious_method_one(double *eps1, double *Z, double *Function, long size, struct com_parameters TR);


