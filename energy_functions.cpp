#include <math.h>
#include "energy_functions.h"
#include <complex>
using namespace std;

void dip_vector(struct parameters *S)
{	/* Функция вычисления прямоугольных координат вектора дипольного момента антенны (дипольное приближение антенны)
		по угловым координатам оси её ДН и типу поляризации.
		
		Входной аргумент - указатель на объект структуры parameters, содержащий параметры передатчика или приёмника
		(в объекте должны быть определены поля структуры S.pz, S.phi0, S.teta0)
		
		Возвращает прямоугольные координаты вектора дипольного момента антенны в системе Oxyz (формулы (2.9))
		в полях структуры S.xp, S.yp, S.zp.
		
		Данную функцию необходимо вызывать при задании полей структуры parameters для приёмника и передатчика
		(перед первым вызовом функций вычисления лучей).
	*/
	
	//вычисление угловых координат (teta_p,phi_p) вектора дипольного момента антенны в системе Oxyz:
	double teta_p, phi_p;
	
	if(S->pz==1) //вектор дипольного момента антенны направлен параллельно поверхности земли
	{
		teta_p=M_PI/2.;
		phi_p=S->phi0+M_PI/2.;
	}
	else //вектор дипольного момента антенны направлен под углом к поверхности земли
	{
		teta_p=S->teta0+M_PI/2.;
		if(teta_p>M_PI || teta_p<0) teta_p=S->teta0-M_PI/2.; //приведение угла teta_p в интервал от 0 до pi
		phi_p=S->phi0;
	}

	//перевод угловых координат вектора дипольного момента антенны в прямоугольные в системе Oxyz (xpt, ypt, zpt) (формулы (2.9)):
	S->xp=sin(teta_p)*cos(phi_p);
	S->yp=sin(teta_p)*sin(phi_p);
	S->zp=cos(teta_p);
	
	return;
}

double ray_delay(double t_nTd,double L, double Td, double *t_,unsigned long *n_,double *dT)
{
	/*вычисление: 
	1)момента времени прихода луча, излученного в момент t_nTd (точного значения t_ и значения в отсчётах n_
	  (точнее, номера ближайшей первой выборки, в которой значение сигнала будет уже не нулевым));
	2)времени распространения луча dt;
	3)интервала dT м/д моментом t_ и следующей за ним ближайшей выборкой
	
	Входные параметры:
	t_nTd - момент излучения луча на n-й итерации
	L - путь луча в текущий момент времени
	Td - период дискретизации сигнала
	*t_ - точное значение момента времени прихода луча
	*n_ - указатель на переменную, в которую будет записана задержка сигнала в отсчётах (относительно нулевого момента времени t=0), т.е. дискретный момент времени прихода луча
	*dT - указатель на переменную, в которую будет записана разность между моментом прихода сигнала и ближайшей за ним выборкой
	
	Выходные параметры:
	dt - возвращаемое значение;
	t_ - записывается по указателю *t_;
	n_ - записывается по указателю *n_;
	dT - записывается по указателю *dT.
	*/

	double dt, flr, dT_T, e=1e-300;
	dt=L/3e+8;//время распространения луча
	*t_=t_nTd+dt;// момент времени прихода луча, излученного в момент t_nTd
	flr=floor((*t_)/Td); //сколько интервалов T целиком укладывается в t_
	*dT=Td*( 1. - ((*t_)/Td-flr) ); //интервал м/д моментом t_ и следующей за ним ближайшей выборкой
	dT_T=(*dT)/Td; //относительное значение интервала м/д моментом t_ и следующей за ним ближайшей выборкой
	if(dT_T <= e || 1.-dT_T <= e) (*dT)=0.; //в силу погрешностей в машинных вычислениях здесь применяется сравнение не с нулём, а с очень малым числом e; dT зануляется, если оно очень близко к 0 или к T, это нужно для обеспечения правильных условий выбора в последующих двух строках, а также, для правильного вычисления величины коррекции начальной фазы зеркально отражённого сигнала, т.к. в случае, когда выборка не приходится точно на момент t0, нужно учесть небольшой набег фазы за время dT
	if((*dT)!=0 || dT_T<=e) *n_=(unsigned long)(flr+1.); //выполняется, если изначально (а не после предыдущей строки) 0<dT<T
	else *n_=(unsigned long)flr; //выполняется, если изначально dT=T, т.е. в t_ содержится целое число интервалов длиной T
	return dt;
}

double pn(double xp,double yp,double zp,double r0x,double r0y,double r0z,double *xE,double *yE,double *zE)
{	/* определение направления вектора E линейно поляризованной волны антенны передатчика в направлении
	   на приёмник (или приёмника на передатчик):
	   
	   Входные параметры:
	   xp, yp, zp - прямоугольные координаты вектора дипольного момента антенны в дипольном приближении;
	   r0x, r0y, r0z - единичный вектор направления (от передатчика на приёмник или от приёмника на передатчик)
	   *xE, *yE, *zE - указатели на переменные, в которые будут записаны координаты вектора E
	   
	   Выходные параметры:
	   длина вектора E - возвращаемое значение
	   xE, yE, zE - координаты вектора E - записываются по указателям *xE, *yE, *zE
	*/

	double pn_=xp*r0x+yp*r0y+zp*r0z; //скалярное произведение (n0,p0) в ф-ле (2.18)
	*xE=xp-pn_*r0x; //← координаты вектора E
	*yE=yp-pn_*r0y; //←
	*zE=zp-pn_*r0z; //←
	return sqrt(pow2(*xE)+pow2(*yE)+pow2(*zE)); //длина вектора E
}


double doppler(char dr,double t_nTd,struct coords C,struct parameters T,struct parameters R,
							 struct com_parameters TR,double *d_2,double *Ld,double *Lr)
{	/*	Функция для вычисления доплеровского сдвига частоты в прямом и зеркально отражённом лучах;
		также вычисляет текущую длину пути этих лучей.
		
		Входные параметры:
		dr - может принимать значения 'd' или 'r' - в зависимости от того, для какого луча нужно вычислить
			 доплеровский сдвиг ('d' - для прямого, 'r' - для зеркально отражённого);
		t_nTd - момент излучения луча на n-й итерации;
		C - объект структуры coords, содержащий интерполированные координаты приёмника и передатчика в текущий момент
			времени t_nTd;
		T - объект структуры parameters, содержащий параметры передатчика;
		R - объект структуры parameters, содержащий параметры приёмника;
		TR - объект структуры com_parameters, содержащий параметры, общие для приёмника и передатчика;
		*d_2 - указатель на переменную, в которую записывается расстояние в проекции на Землю между приёмником
			   и передатчиком в текущий момент времени t_nTd;
	    *Ld - указатель на переменную, в которую записывается расстояние между приёмником и передатчиком
			  (длина прямого луча) в текущий момент времени t_nTd;
		*Lr - указатель на переменную, в которую записывается длина зеркального луча в текущий момент времени t_nTd;
		
		Выходные параметры:
		доплеровский сдвиг частоты луча - возвращаемое значение;
		d_2 - расстояние в проекции на Землю между приёмником и передатчиком в текущий момент времени t_nTd
			  (записывается по указателю *d_2);
	    Ld - расстояние между приёмником и передатчиком (длина прямого луча) в текущий момент времени t_nTd
			  (записывается по указателю *Ld);
	    Lr - длина зеркального луча в текущий момент времени t_nTd (записывается по указателю *Lr);
	*/
	
	double x1t,y1t,z1t,x1r,y1r,z1r,d_2_,L1,Tfc;
	
	unsigned long n_=(unsigned long)ROU(t_nTd/TR.Td);//округлённое целое кол-во периодов дискретизации в интервале времени от 0 до t_nTd
			
	// вычисление доплеровского сдвига частоты для прямого и зеркально отражённого лучей:
	x1t=C.xt+T.Vx[n_]*TR.Td;//координаты передатчика (в СК Oxyz), которые будут через один период дискретизации Td
	y1t=C.yt+T.Vy[n_]*TR.Td;
	z1t=C.zt+T.Vz[n_]*TR.Td;
	
	x1r=C.xr+R.Vx[n_]*TR.Td;//координаты приёмника (в СК Oxyz), которые будут через один период дискретизации Td
	y1r=C.yr+R.Vy[n_]*TR.Td;
	z1r=C.zr+R.Vz[n_]*TR.Td;
	
	d_2_=pow2(x1r-x1t)+pow2(y1r-y1t);
	if(dr=='d') L1=sqrt(d_2_+pow2(z1r-z1t));//расстояние между передатчиком и приёмником, которое будет через один период дискретизации Td
	else L1=sqrt(d_2_+pow2(z1r+z1t));//длина пути зеркального луча между передатчиком и приёмником, которая будет через один период дискретизации Td
	
	*d_2=pow2(C.xr-C.xt)+pow2(C.yr-C.yt);
	*Ld=sqrt(*d_2+pow2(C.zr-C.zt));//текущее расстояние между передатчиком и приёмником
	*Lr=sqrt(*d_2+pow2(C.zr+C.zt));//текущая длина пути зеркального луча между передатчиком и приёмником
	
	//вычисления доплеровских смещений частоты (в Гц) для прямого и зеркального лучей в соостветствии с ф-лой (2.25):
	Tfc=1./TR.Td/TR.lambda;
	if(dr=='d') return (*Ld-L1)*Tfc;
	else return (*Lr-L1)*Tfc;
}

double directional_diagram(double x,double y,double z,double teta0,double phi0,double dteta,
						   double dphi,double n_phi,double n_teta, double C_phi,double C_teta,double C)
/* - функция, вычисляющая нормированное значение ДН по мощности в заданном направлении
Возвращает значение ДН или код ошибки 2, если все координаты вектора направления нулевые

Описание аргументов функции:
x,y,z - координаты вектора направления в СК (Oxyz)", в котором необходимо вычислить величину ДН
		(если координаты даны в СК Oxyz, то нужно вместо x передать в функцию z, а вместо z - x)
teta0, phi0 - угловые координаты оси ДН в сферической СК, связанной либо с передатчиком, либо с приёмником (задаются в радианах)
dteta, dphi - полуширины ДН по соответствующим углам (задаются в радианах)
n_phi, n_teta - показатели степеней, в которые возводятся функции, моделирующие ДН в соответствующих плоскостях
				(влияют на форму ДН)
C_phi, C_teta - параметры a и b в формулах (2.11-2.12) в 1-м отчёте и в формулах, приведённых в начале пункта 1.4 4-го отчёта, определяющие уровень, по которому задаются ширины ДН
C - параметр, определяющий минимально возможное значение пространственной ДН, равное C / ( 1 + C ) (см. формулы, приведённе в начале пункта 1.4 4-го отчёта)
*/
{	double x_, y_, z_, r, phi1, teta1, m, f, a1, a2, a3;
	double sin_teta0, cos_teta0, sin_phi0, cos_phi0, phi_,teta_;

	r=sqrt(x*x+y*y+z*z); //длина вектора направления
	if(r==0) return 2; // error: all the input coordinates are zeros
	x/=r; //нормировка вектора направления
	y/=r;
	z/=r;

	//вычисление вспомогательных (несколько раз используемых ниже) переменных для ускорения вычислений:
	sin_teta0=sin(teta0);
	cos_teta0=cos(teta0);
	sin_phi0=sin(phi0);
	cos_phi0=cos(phi0);
	a2=y*sin_phi0+z*cos_phi0;
	a1=a2*sin_teta0+x*cos_teta0;
	a3=y*cos_phi0-z*sin_phi0;
	
	//перевод координат вектора направления в повёрнутую специальным образом СК по ф-лам (2.16):
	x_=-a2*cos_teta0+x*sin_teta0;
	y_=a3*cos_phi0+a1*sin_phi0;
	z_=-a3*sin_phi0+a1*cos_phi0;
	
	//вычисление в сферической СК, связанной с повёрнутой СК, угловых координат вектора направления по ф-лам (2.8):
	phi_=atan2(y_,z_);
	teta_=acos(x_);
	
	phi1=phi_-phi0; //учёт поворота оси ДН на угол phi0
	teta1=teta_-M_PI/2.; //частичные вычисления по ф-ле (2.11)
	
	//периодизация ф-ии F^2(phi) с периодом 2pi:
	m=floor(fabs(phi1)/2./M_PI);
	if(phi1>2.*M_PI) phi1-=2.*M_PI*m;
	else if(phi1<-2.*M_PI) phi1+=2.*M_PI*m;
	if(phi1<-M_PI) phi1+=2.*M_PI;
	else if(phi1>M_PI) phi1-=2.*M_PI;
	
	//частичные вычисления по ф-ле (2.11):
	if(phi1!=0.)	f=pow(fabs(dphi/C_phi/phi1*sin(C_phi*phi1/dphi)),n_phi);
	else f=1.; //если phi1=0, то будет неопределённость sin(x)/x при x=0, которая равна 1
		
	//частичные вычисления по ф-ле (2.11):
	if(teta1!=0.) f*=pow(fabs(dteta/C_teta/teta1*sin(C_teta*teta1/dteta)*cos(teta1)),n_teta); //если teta1=0, то будет неопределённость sin(x)/x при x=0, которая равна 1
	return (C+f)/(1.+C); //учёт параметра C, и возвращение значения ДН
}

double abs_z3(double *teta,char pz,double eps,double sigma,double f,double sb,double cb)
/* функция вычисляет модуль и аргумент комплексного коэффициента отражения Френеля для TE- или TM-компоненты поля:
   возвращает модуль, а значение аргумента помещает в переменную, на которую указывает указатель *teta
   pz - определяет тип эл/м волны (компоненту поля): pz=0 для ортогональной плоскости падения составляющей эл. поля (TE-компонента), pz=1 - для лежащей в плоскости падения (TМ-компонента)
   eps - диэлектрическая проницаемость среды, от которой происходит отражение
   sigma - проводимость среды
   f - частота эл/м волны
   sb, cb - синус и косинус угла скольжения падающего луча
*/
{
	//структура для описания комплексного числа:
	struct complex
	{
		double x, y;
	};
	
	struct complex a, z, sa, z1, z2, z3;
	double arg_z_2, abs_z_s, g;

	//вычисление вспомогательных переменных:
	a.x=eps;
	a.y=60.*sigma*3e+8/f;
	z.x=a.x-cb*cb;
	z.y=a.y;
	arg_z_2=atan2(z.y, z.x)/2.;
	abs_z_s=sqrt(sqrt(z.x*z.x+z.y*z.y));
	sa.x=abs_z_s*cos(arg_z_2);
	sa.y=abs_z_s*sin(arg_z_2);
	
	//вычисления по ф-лам (1.18):
	if(pz==1)
	{
		z1.x=a.x*sb-sa.x;
		z1.y=a.y*sb-sa.y;
		
		z2.x=a.x*sb+sa.x;
		z2.y=a.y*sb+sa.y;
	}
	else
	{
		z1.x=sb-sa.x;
		z1.y=-sa.y;
		
		z2.x=sb+sa.x;
		z2.y=sa.y;
	}

	g=z2.x*z2.x+z2.y*z2.y;
	z3.x=(z1.x*z2.x+z1.y*z2.y)/g;
	z3.y=(z1.y*z2.x-z1.x*z2.y)/g;

	*teta=atan2(z3.y, z3.x);
	return sqrt(z3.x*z3.x+z3.y*z3.y);
}

int midpoint(double t_nTd,struct parameters *T,struct parameters *R,struct com_parameters *TR,struct coords *C)
{	/*вычисление координат приёмника и передатчика в момент времени t_nTd (момент излучения на n-й итерации)
		Входные параметры:
		t_nTd - момент излучения луча на n-й итерации
		*T - указатель на объект структуры parameters, содержащий параметры передатчика
		*R - указатель на объект структуры parameters, содержащий параметры приёмника
		*TR - указатель на объект структуры com_parameters, содержащий параметры модели, общие для передатчика и приёмника
		*C - указатель на объект структуры coords, в который будут записаны вычисляемые координаты приёмника и
			 передатчика в момент времени t_nTd
		Выходные параметры:
		возвращаемое значение - код ошибки (0 - нет ошибки, 10 - выход за пределы моделируемого интервала времени
			 для заданных n и t)
		C.xt, C.yt, C.zt, C.xr, C.yr, C.zr - координаты приёмника и передатчика в момент времени t_nTd
			 (в объекте C структуры cords)
	*/

	unsigned long n_=(unsigned long)floor(t_nTd/TR->Td);//целое кол-во периодов дискретизации в интервале времени от 0 до t_nTd (момент излучения на n-й итерации)
	if(n_>=TR->N) return 10;//выход за пределы моделируемого интервала
	double dt=t_nTd-(double)n_*TR->Td;
	C->xt = T->x[n_] + T->Vx[n_]*dt;
	C->yt = T->y[n_] + T->Vy[n_]*dt;
	C->zt = T->z[n_] + T->Vz[n_]*dt;
	C->xr = R->x[n_] + R->Vx[n_]*dt;
	C->yr = R->y[n_] + R->Vy[n_]*dt;
	C->zr = R->z[n_] + R->Vz[n_]*dt;
	return 0;
}


long signal(long N, double ampl, double Dfrec, double F_d, double T, double t, double *amplituda, double *phasa, double *frec)
{
	double Cur_n, NumOfPeriod;
	double pi = 3.1415;
	double Period = T*F_d;
	double Period_Impuls = t*F_d;
	unsigned long  int n;
	complex <double> One(0.0, 1.0), omega(2*pi*10, 0);


	for(n = 0; n < N; n++)
	{
		NumOfPeriod = floor(n/(Period));
		Cur_n = n - NumOfPeriod*Period;

		if(Cur_n<=Period_Impuls)
		{
			*(phasa + n) = ((Dfrec*Cur_n*Cur_n)/(4*pi*Period_Impuls) - Dfrec*Cur_n/(4*pi))/F_d;
			*(amplituda + n) = ampl;
		}
		else
		{
			*(phasa + n) = 0;
			*(amplituda + n) = 0;
		}

	
		complex <double> phi(phasa[n],0);
		complex <double> a(amplituda[n],0);
		complex <double> tt(n/F_d,0);
		frec[n] = arg(a*exp(omega*One*tt)*exp(phi*One));

	}

	return 0;
}

// двумерная функции для вычисления направления на цель методом максимального правдоподобия
int max_specious_method(double *eps1, double *eps2, double *Z, double *Function, long size, struct com_parameters TR){
	
	int dim = 100;
	double phi_temp2 = 0, phi_temp1 = 0, phi_min1  = 0, phi_min2  = 0;
	double Function_temp[100*100] = {0};
	complex <long double> C[2][2], inv_C[2][2], M[2][2], Y1[5][1], Y2[1][5], One(0, 1);
	complex <long double> S1[5], S2[5], F1[5][2], F2[2][5], FF[2][2], inv_FF[2][2];
	complex <long double> M1[5][2], M2[5][5], M3[1][5], Functional(0.0,0.0);
	int K = 0, L = 0;

	double b = -1e+100, f_max = -1e+100;

	for (int nn = 0; nn<5; nn++)
	{
		Y1[nn][0] = (complex <long double>)Z[2*nn] + One*(complex <long double>)Z[2*nn + 1];
		Y2[0][nn] = conj((complex <long double>)Z[2*nn] + One*(complex <long double>)Z[2*nn + 1]);	
	}

	for (phi_temp1 = 0.; phi_temp1<(M_PI/2); phi_temp1 += (M_PI/2)/(dim-1))
	{
		K = 0;

		int t = 0, t1 = 0, t2 = 0, t3 = 0;
		for(int nn = 0; nn<5; nn++){
			S1[nn] = exp((complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp1)*(complex <long double>)nn);
		}				

		for (phi_temp2 = -M_PI/2; phi_temp2<0; phi_temp2 += (M_PI/2)/(dim-1))
		{
				for(int nn = 0; nn<5; nn++){
					S2[nn] = exp((complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp2)*(complex <long double>)nn);
				}

			for(int nn = 0; nn<5; nn++)
				{
					F1[nn][0] = S1[nn];
					F1[nn][1] = S2[nn];
					F2[0][nn] = conj(S1[nn]);
					F2[1][nn] = conj(S2[nn]);
				}
			for (int i = 0; i<2; i++){
				for (int j = 0; j<2; j++){
						FF[i][j] = (complex <long double>)0 + (complex <long double>)0*One;
						for(int r = 0; r<5; r++){
							FF[i][j] += F2[i][r]*F1[r][j];
					}
				}
			}

			if ((FF[0][0]*FF[1][1] - FF[0][1]*FF[1][0])!= (complex <long double>)(0,0))
			{
				inv_FF[0][0] = FF[1][1];
				inv_FF[1][0] = -FF[1][0];
				inv_FF[0][1] = -FF[0][1];
				inv_FF[1][1] = FF[0][0];

				for (int i = 0; i<2; i++){
					for (int j = 0; j<2; j++){
							inv_FF[i][j] /= (FF[0][0]*FF[1][1] - FF[0][1]*FF[1][0]);
						}
				}

			} else {
				for (int i = 0; i<2; i++){
					for (int j = 0; j<2; j++){
						inv_FF[i][j] = (complex <long double>)0;
					}
				}
			}

			for (int i = 0; i<5; i++){
				for (int j = 0; j<2; j++){
						M1[i][j] = (complex <long double>)0 + (complex <long double>)0*One;
						for(int r = 0; r<2; r++){
							M1[i][j] += F1[i][r]*inv_FF[r][j];
						}
				}
			}


			for (int i = 0; i<5; i++){
				for (int j = 0; j<5; j++){
						M2[i][j] = (complex <long double>)0 + (complex <long double>)0*One;
						for(int r = 0; r<2; r++){
							M2[i][j] += M1[i][r]*F2[r][j];
						}
				}
			}

			for (int i = 0; i<1; i++){
				for (int j = 0; j<5; j++){
						M3[i][j] = (complex <long double>)0 + (complex <long double>)0*One;
						for(int r = 0; r<5; r++){
							M3[i][j] += Y2[i][r]*M2[r][j];
						}
				}
			}

			for (int i = 0; i<1; i++){
				for (int j = 0; j<1; j++){
						Functional = (complex <long double>)0 + (complex <long double>)0*One;
						for(int r = 0; r<5; r++){
							Functional += M3[i][r]*Y1[r][j];
						}
				}
			}

			Function_temp[K + L*(dim)] = abs(Functional);
			
			K++;
		}

		L++;
	}
	
	K = 0;
	L = 0;

	for ( phi_temp1 = 0.; phi_temp1<(M_PI/2); phi_temp1 += (M_PI/2)/(dim-1)){
		K = 0;
		for (phi_temp2 = -M_PI/2; phi_temp2<(0); phi_temp2 += (M_PI/2)/(dim-1)){
			if ((K!=0)&&((L!=0))&&(K<(dim-1))&&(L<(dim-1)))
			{
				if((Function_temp[K + L*(dim)]>Function_temp[K + L*(dim)-1])&&(Function_temp[K + L*(dim)]>Function_temp[K + L*(dim)+1])&&(Function_temp[K + L*(dim)]>Function_temp[K + (L-1)*(dim)])&&(Function_temp[K + L*(dim)]>Function_temp[K + (L+1)*(dim)]))
				{
					if (f_max < Function_temp[K + L*(dim)]){
					
						f_max = Function_temp[K + L*(dim)];
						phi_min1 = phi_temp1;
						phi_min2 = phi_temp2;		
						if (phi_min1<phi_min2){
						phi_min1 = phi_temp2;
						phi_min2 = phi_temp1;
						}

					}
				
				}
			}

			Function[(int)floor((double)(K*size)/dim) + (int)floor((double)(L*size)/dim)*size] = Function_temp[K + L*(dim)];			
			K++;
		}
		L++;
	}



	if (abs(phi_min2 - phi_min1) < 2*M_PI/(dim-1))
	{
		phi_min1 = (phi_min1 + phi_min2)/2;
		phi_min2 = phi_min1;
	
	}

return 0;
}

// одномерная функции для вычисления направления на цель методом максимального правдоподобия
int max_specious_method_one(double *eps, double *Z, double *Function, long size, struct com_parameters TR)
{
	int dim = 500;
	double Function_temp[500];
	double phi_temp = 0, phi_min  = 0;
	double a = 0;
	complex <long double> AA = (0,0), BB = (0,0), CC = (0,0), One(0,1);
	long double Z2_min = -1e+100;
	double S[10], SS[10];
	
	long K = 0;
	int t = 0;

	for (phi_temp = 0; phi_temp < M_PI/2; phi_temp += M_PI/2/(dim-1))
	{
		AA = (0,0); BB = (0,0); CC = (0,0);
		t = 0;
		for(int nn = 0; nn<5; nn++)
			{
				if ((2*M_PI*TR.d2/(TR.lambda)*sin(phi_temp)*nn) > (2*M_PI)){
					t = 1;	
				}
				else {	
					S[2*nn] = real(exp((complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn));
					S[2*nn + 1] = imag(exp((complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn));
					SS[2*nn] = real(conj(exp((complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn)));
					SS[2*nn + 1] = imag(conj(exp((complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn)));

				}
		}

		if (t == 1){	
			for(int nn = 0; nn<5; nn++)
				{
				S[2*nn] = real(exp((complex <long double>)(0*M_PI) + (complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn));
				S[2*nn + 1] = imag(exp((complex <long double>)(0*M_PI) + (complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn));
				SS[2*nn] = real(conj(exp((complex <long double>)(0*M_PI) + (complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn)));
				SS[2*nn + 1] = imag(conj(exp((complex <long double>)(0*M_PI) + (complex <long double>)(-2)*(complex <long double>)M_PI*One*((complex <long double>)TR.d2/(complex <long double>)(TR.lambda))*sin((complex <long double>)phi_temp)*(complex <long double>)nn)));
			}
		}


		for (int nn = 0; nn<5; nn++)
		{	
			AA += ((complex <long double>)Z[2*nn] + One*(complex <long double>)Z[2*nn + 1])*((complex <long double>)Z[2*nn] + One*(complex <long double>)Z[2*nn + 1]);
			BB += ((complex <long double>)Z[2*nn] + One*(complex <long double>)Z[2*nn + 1])*((complex <long double>)SS[2*nn] + One*(complex <long double>)SS[2*nn + 1]);
			CC += ((complex <long double>)S[2*nn] + One*(complex <long double>)S[2*nn + 1])*((complex <long double>)SS[2*nn] + One*(complex <long double>)SS[2*nn + 1]);
		}

		Function_temp[K] = norm(BB/sqrt(CC));
		K++;
	}
	K = 0;
	phi_temp = 0;
	for (phi_temp = 0; phi_temp < M_PI/2; phi_temp += M_PI/2/(dim-1))
	
	{
		if (K>0 && K < (dim-1))
		{
			if ((Function_temp[K] > Function_temp[K-1])&&(Function_temp[K] > Function_temp[K+1]))
			{
				if (Z2_min < Function_temp[K])
				{
					Z2_min = Function_temp[K];
					phi_min = phi_temp;
				}
			}
		}
	
		Function[(int)floor((double)(K*size)/dim)] = Function_temp[K];	
		K++;
	}

	*eps  = phi_min*180/M_PI;

return 0;
}

long DN(double phi, double* ves, long size_ves ,double  lamda, double  d1, double* azimyt,
				double eps, double* amp, double* faza, long size_amp,  double  d2, double* ygol_mesta)
{	
	complex <double> Sum_phi(0.0, 0.0), Sum2_phi(0.0, 0.0), Sum_eps(0.0, 0.0), Sum2_eps(0.0, 0.0);
	complex <double> One(0.0, 1.0), One2(0.0, -1.0);
	complex <double> deg(M_PI/180, 0.0);
		
	for(int i = 0; i < size_ves; i++) 	
		Sum_phi += ves[i]*exp(-2*M_PI*(d1/lamda)*One*((double)i-((double)size_ves-1)/2)*sin(phi*M_PI/180));

	for(int n = 0; n < size_amp; n++)
		Sum_eps += exp(One*faza[n]*deg)*amp[n]*exp(2*M_PI*(d2/lamda)*One*((double)n-((double)size_amp-1)/2)*sin(eps*M_PI/180));

	Sum_phi *= cos((d1/lamda)*phi*M_PI/180)*cos((d1/lamda)*phi*M_PI/180);
	Sum_eps *= cos(eps*M_PI/180);

	azimyt[0] = abs(Sum_phi);
	azimyt[1] = arg(Sum_phi);
	
	ygol_mesta[0] = abs(Sum_eps);
	ygol_mesta[1] = arg(Sum_eps);

return 0;
}

int direct_ray_power(double *pfd,double *P,double *phi,double *phi_dT,double *t_,unsigned long *n_,double *F2t,
					 double *F2r,double *K_d,double *dFdop_d,double *dt,struct parameters T,struct parameters R,
					 struct com_parameters TR,double t_nTd,double phi0,double P0, int chanel, double *mas, double *D, double *debug)
{
	/*
		Расчитывает текущую точку прямого луча в дискретный момент времени t_nTd (мощность, плотность потока мощности,
		начальную фазу, доплеровский сдвиг частоты,время распространения луча, коэффициент рассогласования из-за несовпадения
		поляризаций, значения диаграмм направленности приёмной и передающей антенн в направлении прямого луча, точный
		и дискретный моменты прихода луча,набег начальной фазы между этими моментами).

		Входные параметры:
		*pfd - указатель на переменную, в которую будет записано значение плотности потока мощности прямого луча
			в точке приёма в момент прихода t_, излученного в момент t_nTd;
		*P - указатель на переменную, в которую будет записано значение мощности прямого луча в точке приёма в
			момент прихода t_, излученного в момент t_nTd;
		*phi - указатель на переменную, в которую будет записано значение начальной фазы прямого луча в точке приёма в
			момент прихода прямого луча t_, излученного в момент t_nTd; начальная фаза вычисляется с учётом влияния
			доплеровского сдвига частоты;
		*phi_dT - указатель на переменную, в которую будет записано значение набега начальной фазы между точным t_ и 
			дискретным n_ (отображаемом на графиках) моментами прихода прямого луча в точку приёма;
		*t_ - указатель на переменную, в которую будет записано значение точного момента прихода прямого луча
			в точку приёма;
		*n_ - указатель на переменную, в которую будет записано значение дискретного (ближайшего, следующего за точным)
			момента прихода прямого луча в точку приёма;
		*F2t - указатель на переменную, в которую будет записано значение диаграммы направленности антенны передатчика
			по мощности в направлении прямого луча в момент излучения t_nTd;
		*F2r - указатель на переменную, в которую будет записано значение диаграммы направленности антенны приёмника
			по мощности в направлении прямого луча в точный момент его прихода t_;
		*K_d - указатель на переменную, в которую будет записано значение поляризационного коэффициента рассогласования
			по амплитуде для прямого луча в точный момент его прихода t_;
		*dFdop_d - указатель на переменную, в которую будет записано значение доплеровского сдвига частоты прямого луча
			в точный момент его прихода t_;
		*dt - указатель на переменную, в которую будет записано значение времени распространения прямого луча,
			излученного в момент t_nTd;
		T - объект структуры parameters, содержащий параметры передатчика;
		R - объект структуры parameters, содержащий параметры приёмника;
		TR - объект структуры com_parameters, содержащий параметры модели, общие для приёмника и передатчика;
		t_nTd - точный момент излучения прямого луча на n-й итерации;
		phi0 - начальная фаза излучаемого передатчиком сигнала в точный момент начала его излучения t;
		P0 - мощность, подводимая к антенне передатчика;
				
		Выходные параметры:
		возвращаемое значение - код ошибки (0 - нет ошибки, 2 - нулевое расстояние м/д передатчиком и приёмником,
			10 - выход за пределы моделируемого интервала времени для заданных n и t);
		все переменные, на которые были переданы указатели в качестве входных параметров (см. описание входных параметров).
	*/

	double r0x,r0y,r0z,r,dd,r1,r2,xEt,yEt,zEt,xEr,yEr,zEr,d_2,Ld,Lr,dT,err_no;
	double phasa1 = 0, phasa2 = 0;	
	struct coords C;
	
	err_no= midpoint(t_nTd,&T,&R,&TR,&C);

	if(err_no==10) return 10;// выход сигнала за пределы моделируемого интервала
	
	*pfd=P0*T.G/4./M_PI;

	r0x=C.xr-C.xt; //вычисление вектора r0, направленного от передатчика к приёмнику, в СК Oxyz
	r0y=C.yr-C.yt;
	r0z=C.zr-C.zt;
	r2=pow2(r0x)+pow2(r0y)+pow2(r0z);// квадрат длины вектора r0
	if(r2==0) return 2; //error: r0=0 - нулевое расстояние м/д передатчиком и приёмником
	(*pfd)/=r2; //  /r0^2 - частичное вычисление по ф-лам (2.23) и (2.24)
	//if(T.i==0) dd = My_directional_diagram(TR,r0x,r0y,r0z, T.teta0,T.phi0, T.T_ves, T.size_ves, TR.d1, T.T_amp, T.T_faza, T.T_size_amp, TR.d2, &phasa1, T.MaxDN_T, chanel, mas); //вычисление ДН антенны передатчика в направлении приёмника, если антенна передатчика не изотропная
	
	if(T.i==0)
	{
		dd=My_directional_diagram_3(r0x,r0y,r0z, T.teta0,T.phi0,T.dteta,T.dphi,T.n_phi,T.n_teta,T.C_phi, T.C_teta, T.C, debug);
		dd = sqrt(dd);
//		dd = 1.;
	}
	else dd=1.;// если антенна передатчика изотропная, то ДН=1


	if(dd==2) return 2;// невозможно вычислить ДН, если длина вектора r0 нулевая
	(*pfd)*=(dd*dd); // частичное вычисление по ф-лам (2.23) и (2.24)
	*F2t=(dd*dd); // передача значения нормированной ДН антенны передатчика по мощности в направлении прямого луча
	*dFdop_d=doppler('d',t_nTd,C,T,R,TR,&d_2,&Ld,&Lr);

	*P=(*pfd)*R.G*pow2(TR.lambda)/4./M_PI; //умножение плотности потока мощности на эффективную площадь приёмной антенны
	if(R.i==0){
		dd = directional_diagram(-r0x,-r0y,-r0z, R.teta0,R.phi0,R.dteta,R.dphi,R.n_phi,R.n_teta,R.C_phi, R.C_teta, R.C, debug);
		dd = sqrt(dd);
		directional_diagram(TR, -r0x,-r0y,-r0z, R.teta0,R.phi0, R.R_ves, R.size_ves, TR.d1, R.R_amp,
			R.R_faza, R.R_size_amp, TR.d2, &phasa2, R.MaxDN_R, chanel, mas, debug); // вычисление ДН антенны приёмника в направлении передатчика, если антенна приёмника не изотропная
																					//	d=1.;
	}
	else dd=1.;
	if(dd==2) return 2;// невозможно вычислить ДН, если длина вектора r0 нулевая
	(*P)*=(dd*dd);// частичное вычисление по ф-лам (2.23) и (2.24)
	*D = (dd*dd);
	*F2r=(dd*dd); // передача значения нормированной ДН антенны приёмника по мощности в направлении прямого луча
	
	r=sqrt(r2);
	r0x/=r;// нормировка вектора r0
	r0y/=r;
	r0z/=r;

	// определение направления вектора E линейно поляризованной волны антенны передатчика в направлении на приёмник:
	r1=pn(T.xp,T.yp,T.zp,r0x,r0y,r0z,&xEt,&yEt,&zEt); //длина вектора E
	
	// вычисление плоскости поляризации антенны приёмника (в системе Oxyz в направлении на передатчик):
	r2=pn(R.xp,R.yp,R.zp,-r0x,-r0y,-r0z,&xEr,&yEr,&zEr); //длина вектора поляризации
	
	//нормировка вектороа E линейно поляризованной волны антенны передатчика и вектора поляризации антенны приёмника:
	if(r1!=0 && r2!=0)
	{
		xEt/=r1;
		yEt/=r1;
		zEt/=r1;
		
		xEr/=r2;
		yEr/=r2;
		zEr/=r2;
		
		*K_d=fabs(xEt*xEr+yEt*yEr+zEt*zEr); //вычисление амплитудного коэффициента рассогласования по поляризации (2.19)
	}

	(*P)*=pow2(*K_d); //учёт поляризационного рассогласования

	phi0  += (phasa1 + phasa2);
	*dt=ray_delay(t_nTd, Ld, TR.Td, t_, n_, &dT);//вычисление момента прихода прямого сигнала (в секундах и в отсчётах),
	//а также времени распротсранения луча dt и промежутка времени dT между моментом прихода луча и ближайшей за ним выборкой
	if(*n_>=TR.N) return 10;//выход за пределы моделируемого интервала		
	*phi_dT=2.*M_PI*(*dFdop_d)*dT;//набег начальной фазы за промежуток времени dT между моментом прихода луча и ближайшей за ним выборкой
	*phi=phi0-2*M_PI/TR.lambda*Ld;//начальная фаза прямого сигнала в точке приёма в момент времени t_ (момент прихода луча в точку приёма, излучённого в момент t+n*TR.Td)

	return 0;
}


int reflected_ray_power(double *pfd1,double *pfd2,double *teta1,double *teta2,double *P,double *phi,double *phi_,
						double *phi_dT,double *t_,unsigned long *n_,double *F2t,double *F2r,double *K_r,
						double *dFdop_r,double *dt,struct parameters T,struct parameters R,struct com_parameters TR,
						double t_nTd,double phi0,double P0, int chanel, double *mas, double *D, double *debug)
{
	/*
		Расчитывает текущую точку квазизеркально отражённого луча в дискретный момент времени t_nTd (мощность,
		плотность потока мощности для TE- и TM-компонент, сдвиги фаз при отражении от подстилающей поверхности
		для TE- и TM-компонент, начальную фазу с учётом и без учёта сдвига при отражении, доплеровский сдвиг частоты,
		время распространения луча, коэффициент рассогласования из-за несовпадения поляризаций, значения диаграмм
		направленности приёмной и передающей антенн в направлении прямого луча, точный и дискретный моменты прихода
		луча, набег начальной фазы между этими моментами).

		Входные параметры:
		*pfd1 - указатель на переменную, в которую будет записано значение плотности потока мощности TE-компоненты
			(горизонтальная поляризация) квазизеркально отражённого луча в точке приёма в момент прихода t_,
			излученного в момент t_nTd;
		*pfd2 - указатель на переменную, в которую будет записано значение плотности потока мощности TM-компоненты
			(вертикальная поляризация) квазизеркально отражённого луча в точке приёма в момент прихода t_,
			излученного в момент t_nTd;
		*teta1 - указатель на переменную, в которую будет записано значение сдвига фазы TE-компоненты квазизеркально
			отражённого луча, обусловленного отражением от подстилающей поверхности;
		*teta2 - указатель на переменную, в которую будет записано значение сдвига фазы TM-компоненты квазизеркально
			отражённого луча, обусловленного отражением от подстилающей поверхности;
		*P - указатель на переменную, в которую будет записано значение мощности зеркального луча в точке приёма в
			момент прихода зеркального луча t_, излученного в момент t_nTd;
		*phi - указатель на переменную, в которую будет записано значение начальной фазы зеркального луча в точке
			приёма в момент прихода зеркального луча t_, излученного в момент t_nTd с учётом сдвига фазы при
			отражении от подстилающей поверхности и с учётом доплеровского сдвига частоты;
		*phi_ - указатель на переменную, в которую будет записано значение начальной фазы зеркального луча в точке
			приёма в момент прихода зеркального луча t_, излученного в момент t_nTd без учёта сдвига фазы при
			отражении от подстилающей поверхности, но с учётом доплеровского сдвига	частоты;
		*phi_dT - указатель на переменную, в которую будет записано значение набега начальной фазы между точным t_ и 
			дискретным n_ (отображаемом на графиках) моментами прихода зеркального луча в точку приёма;
		*t_ - указатель на переменную, в которую будет записано значение точного момента прихода зеркального луча
			в точку приёма;
		*n_ - указатель на переменную, в которую будет записано значение дискретного (ближайшего, следующего за точным)
			момента прихода зеркального луча в точку приёма;
		*F2t - указатель на переменную, в которую будет записано значение диаграммы направленности антенны передатчика
			по мощности в направлении зеркального луча в момент излучения t_nTd;
		*F2r - указатель на переменную, в которую будет записано значение диаграммы направленности антенны приёмника
			по мощности в направлении прямого зеркального в точный момент его прихода t_;
		*K_r - указатель на переменную, в которую будет записано значение поляризационного коэффициента рассогласования
			по амплитуде для прямого луча в точный момент его прихода t_;
		*dFdop_r - указатель на переменную, в которую будет записано значение доплеровского сдвига частоты зеркального
			луча в точный момент его прихода t_;
		*dt - указатель на переменную, в которую будет записано значение времени распространения зеркального луча,
			излученного в момент t_nTd;
		T - объект структуры parameters, содержащий параметры передатчика;
		R - объект структуры parameters, содержащий параметры приёмника;
		TR - объект структуры com_parameters, содержащий параметры модели, общие для приёмника и передатчика;
		t_nTd - точный момент излучения зеркального луча на n-й итерации;
		phi0 - начальная фаза излучаемого передатчиком сигнала в точный момент начала его излучения t;
		P0 - мощность, подводимая к антенне передатчика;
		n - номер текущей итерации алгоритма моделирования
		
		Выходные параметры:
		возвращаемое значение - код ошибки (0 - нет ошибки, 2 - нулевое расстояние м/д передатчиком и приёмником,
			10 - выход за пределы моделируемого интервала времени для заданного t_nTd);
		все переменные, на которые были переданы указатели в качестве входных параметров (см. описание входных параметров).
	*/

	double d_2,Ld,Lr,x,la,x0,y0,r1x,r1y,r1z,r1,r2x,r2y,r2z,r2,r,dd,nExt,nEyt,nEzt,nExr,nEyr,nEzr;
	double n0x,n0y,E1,E2,sb,cb,cos_alfa,sin_alfa,abs_gamma1,abs_gamma2,dphi,dT;
	int err_no;
	double phasa2 = 0;
	
	struct coords C;
	
	err_no=midpoint(t_nTd,&T,&R,&TR,&C);//вычисление координат приёмника и передатчика в момент t_nTd (момент излучения на n-й итерации)
	if(err_no==10) return 10;// выход сигнала за пределы моделируемого интервала

	*P=P0*T.G/4./M_PI;
				
	*dFdop_r=doppler('r',t_nTd,C,T,R,TR,&d_2,&Ld,&Lr);

	//вычисление координат (x0,y0) точки отражения зеркального луча по формулам на стр. 36 (внизу):
	x=C.zt*sqrt(d_2)/(C.zt+C.zr);
	la=x/(sqrt(d_2)-x);
	x0=(C.xt+la*C.xr)/(1+la);// coordinates of
	y0=(C.yt+la*C.yr)/(1+la);// reflection point
	
	//вычисление вектора r1, направленного от передатчика к точке зеркального отражения, и его длины:
	r1x=x0-C.xt;
	r1y=y0-C.yt;
	r1z=-C.zt;
	r1=sqrt(pow2(r1x)+pow2(r1y)+pow2(r1z));
	
	//вычисление вектора r2, направленного от точки зеркального отражения к приёмника, и его длины:
	r2x=C.xr-x0;
	r2y=C.yr-y0;
	r2z=C.zr;
	r2=sqrt(pow2(r2x)+pow2(r2y)+pow2(r2z));
	
	r=r1+r2; // длина зеркального луча
	
	(*P)/=pow2(r);//частичный расчёт мощности сигнала квазизеркально отражённого луча
	
	if(T.i==0) {
		dd = My_directional_diagram_3(r1x,r1y,r1z,T.teta0,T.phi0,T.dteta,T.dphi,T.n_phi,T.n_teta,T.C_phi,T.C_teta,T.C, debug);
		dd = sqrt(dd);
		//dd=1.;
	}
	else dd=1.;// если антенна передатчика изотропная, то ДН=1
	if(dd==2) return 2; // невозможно вычислить ДН, если длина вектора r1 нулевая
	(*P)*=(dd*dd); //частичный расчёт мощности сигнала квазизеркально отражённого луча
	*F2t=(dd*dd); //передача значения нормированной ДН антенны передатчика по мощности в направлении зеркального луча
	
	r1x/=r1;// нормировка вектора r1
	r1y/=r1;
	r1z/=r1;

	// определение направления вектора E линейно поляризованной волны антенны передатчика в направлении на точку зеркального отражения:
	r=pn(T.xp,T.yp,T.zp,r1x,r1y,r1z,&nExt,&nEyt,&nEzt);
								
	double k0; //вспомогательная переменная
	if(r!=0)
	{
		nExt/=r; //нормировка вектора E
		nEyt/=r;
		nEzt/=r;
		k0=1.;
	}
	else k0=0.; //если длина вектора E нулевая, то k0=0 и коэффициент Френеля при отражении не вычисляется 
	
	n0x=r1y;// вычисление вектора нормали к плоскости падения зеркального луча
	n0y=-r1x;
	r=sqrt(pow2(n0x)+pow2(n0y));
	if(r!=0) //если длина вектора нормали ненулевая:
	{
		n0x/=r; //нормировка вектора нормали
		n0y/=r;
		E1=nExt*n0x+nEyt*n0y; //вычисление составляющей вектора E, в падающей волне, ортогональной плоскости падения
	}
	else E1=1.; //если длина вектора нормали ненулевая (очень редкий случай, в этом случае нет плоскости падения, т.к. передатчик и приёмник расположены друг над другом или в одной и той же точке), то для определённости считаем, что E1=1
	
	E2=sqrt(1.-pow2(E1));// вычисление составляющей вектора E, в падающей волне, лежащей в плоскости падения
	if(nEzt<0) E2=-E2;
	
	r2x/=r2; //нормировка вектора r2
	r2y/=r2;
	r2z/=r2;
	
	// вычисление плоскости поляризации антенны приёмника (в системе Oxyz в направлении на точку зеркального отражения):
	r=pn(R.xp,R.yp,R.zp,-r2x,-r2y,-r2z,&nExr,&nEyr,&nEzr);
								
	if(r!=0) //нормировка вектора поляризации, если он ненулевой
	{
		nExr/=r;
		nEyr/=r;
		nEzr/=r;
		k0=1.;
	}
	else k0=0.; //если длина вектора поляризации нулевая, то k0=0 и коэффициент Френеля при отражении не вычисляется 
	
	sb=C.zt/r1; //sin угла скольжения падающего луча

	double em, Gamma; //вспомогательные переменные
	em=4.*M_PI*TR.ds*sb/TR.lambda;
	Gamma=exp(-pow2(em)); //коэффициент ослабления зеркального сигнала из-за неровностей отражающей поверхности
	(*P)*=Gamma; //умножение плотности потока мощности квазизеркально отражённого луча на этот коэф-т
	
	*pfd1=*P;
	*pfd2=*P;
	
	if(k0!=0.) //если плоскость падения существует, то вычисляется комплексный коэффициент отражения Френеля:
	{
		cb=sqrt(1.-pow2(sb)); // cos угла скольжения падающего луча
		abs_gamma1=abs_z3(teta1,0,TR.eps,TR.sigma,TR.f,sb,cb); // вычисление модуля и аргумента (teta1) коэффициента отражения Френеля для ортогональной плоскости падения составляющей электрического поля
		abs_gamma2=abs_z3(teta2,1,TR.eps,TR.sigma,TR.f,sb,cb); // вычисление модуля и аргумента (teta2) коэффициента отражения Френеля для лежащей в плоскости падения составляющей электрического поля

		(*pfd1)*=pow2(E1*abs_gamma1);//плотность потока мощности зеркального луча в точке приёма, приходящаяся на TE-компоненту (горизонтальную поляризацию)
		(*pfd2)*=pow2(E2*abs_gamma2);//плотность потока мощности зеркального луча в точке приёма, приходящаяся на TM-компоненту (вертикальную поляризацию)

		cos_alfa=nExr*n0x+nEyr*n0y; //cos угла альфа (рис. 2.7) между нормалью к плоскости падения и вектором поляризации антенны приёмника
		sin_alfa=sqrt(1.-pow2(cos_alfa)); //sin угла альфа (рис. 2.7)
		if(nEzr<0.) sin_alfa=-sin_alfa;  //sin угла альфа д.б. отрицательным, если z-компонента вектора поляризации <0
			
		dphi=atan2(E1*cos_alfa*abs_gamma1*sin(*teta1)+E2*sin_alfa*abs_gamma2*sin(*teta2),E1*cos_alfa*abs_gamma1*cos(*teta1)+E2*sin_alfa*abs_gamma2*cos(*teta2)); //разность фаз м/д падающим на подстилающую поерхность полем и проекцией отражённого поля на плоскость поляризации приёмной антенны
														
		*K_r=sqrt(pow2(E1*abs_gamma1*cos_alfa)+pow2(E2*abs_gamma2*sin_alfa)+
			2*E1*E2*abs_gamma1*abs_gamma2*sin_alfa*cos_alfa*cos(*teta1-*teta2));// вычисление поляризационного коэффициента рассогласования по амплитуде с учётом коэффициентов Френеля по ф-ле (2.21)
	}
	else
	{
		*teta1=0;// если плоскость падения не существует, то сдвиги фаз TE- и TM-компонент зеркального луча, обусловленные отражением, полагаются нулевыми
		*teta2=0;
	}
				
	*dt=ray_delay(t_nTd, Lr, TR.Td, t_, n_, &dT);//вычисление момента прихода зеркального луча (в секундах и в отсчётах),
	//а также времени распротсранения луча dt и промежутка времени dT между моментом прихода луча и ближайшей за ним выборкой
	if(*n_>=TR.N) return 10;//выход за пределы моделируемого интервала		

	if(R.i==0){
		dd=My_directional_diagram_3(-r2x,-r2y,-r2z,R.teta0,R.phi0,R.dteta,R.dphi,R.n_phi,R.n_teta,R.C_phi,R.C_teta,R.C, debug); //вычисление ДН антенны приёмника в направлении на точку зеркального отражения, если антенна приёмника не изотропная
		My_directional_diagram_2(TR,-r2x,-r2y,-r2z, R.teta0,R.phi0, R.R_ves, R.size_ves, TR.d1, R.R_amp,
			R.R_faza, R.R_size_amp, TR.d2, &phasa2, R.MaxDN_R, chanel, mas, debug);
		dd = sqrt(dd);
//		dd=1.;
	}
	else dd=1.;// если антенна приёмника изотропная, то ДН=1
	if(dd==2) return 2; // невозможно вычислить ДН, если длина вектора r2 нулевая
	(*P)*=(dd*dd); //частичный расчёт мощности сигнала квазизеркально отражённого луча
	*F2r=(dd*dd); //передача значения нормированной ДН антенны приёмника по мощности в направлении зеркального луча
	*D = (dd*dd);			
	(*P)*=pow2(*K_r); //умножение амплитуды зеркального сигнала на поляризационный коэффициент рассогласования
	//если же k0=0, то в качестве K_r берётся его значение в предыдущей точке, вычисленное при предыдущем вызове DLL-модуля из LabVIEW и переданное в DLL-модуль при текущем его вызове из LabVIEW
	
	(*P)*=R.G*pow2(TR.lambda)/4./M_PI; //умножение на эффективную площадь приёмной антенны
	*phi_dT=2.*M_PI*(*dFdop_r)*dT;//набег начальной фазы за промежуток времени dT между моментом прихода луча и ближайшей за ним выборкой
	
	*phi_=phi0-2*M_PI/TR.lambda*Lr;//начальная фаза квазизеркально отражённого сигнала в точке приёма в момент
					// времени t_ (момент прихода луча в точку приёма, излучённого в момент t+n*TR.Td)
					// без учёта сдвига фазы луча при отражении от подстилающей поверхности dphi
	*phi=*phi_+dphi;// начальная фаза квазизеркально отражённого сигнала в точке приёма в момент
					// времени t_ (момент прихода луча в точку приёма, излучённого в момент t+n*TR.Td)
					// с учётом сдвига фазы луча при отражении от подстилающей поверхности dphi
	return 0;
}


int disp_ray_power(int s,double *P_disp,double *sin_mean_,long *N_delta_2,struct parameters T,
				   struct parameters R,struct com_parameters TR,double t_nTd,double P0,unsigned long n,int chanel, double *mas, double *D, double *debug)
{
	
	/*
		Вычисляет мощность рассеянного сигнала для одного момента моделирования n, либо для всех моментов за один вызов,
		если мощность с течением времени не меняется (приёмник и передатчик неподвижны); также вычисляет зависимость
		от времени среднего по всем рассеивающим площадкам значения выражения sin(beta1)+sin(beta2) и кол-во
		рассеивающих площадок в текущий момент времени моделирования n

		Входные параметры:
		s - параметр, указывающий на режим работы функции: при s=0 мощность рассеянного сигнала вычисляется для
			одного момента моделирования n (для каждого n необходимо делать отдельный вызов); данный режим используется
			при подвижных приёмнике и (или) передатчике, когда мощность рассеянного сигнала будет зависеть от времени;
			при s=1 мощность рассеянного сигнала вычисляется для всех моментов времени моделирования n за один вызов
			(при n=0); данный режим используется, если мощность с течением времени не меняется (приёмник и передатчик
			неподвижны);
		*P_disp - указатель на массив значений мощности рассеянного сигнала;
		*sin_mean_ - указатель на массив значений среднего по всем рассеивающим площадкам значения выражения
			sin(beta1)+sin(beta2);
		*N_delta_2 - указатель на переменную, в которую записывается кол-во рассеивающих площадок в текущий момент
			времени моделирования n (для сигнала, излученного в момент времени t+n*Td);
		T - объект структуры parameters, содержащий параметры передатчика;
		R - объект структуры parameters, содержащий параметры приёмника;
		TR - объект структуры com_parameters, содержащий параметры модели, общие для приёмника и передатчика;
		t_nTd - точный момент излучения лучей, образующих рассеянный сигнал на n-й итерации;
		P0 - мощность, подводимая к антенне передатчика;
		
		Выходные параметры:
		возвращаемое значение - код ошибки (0 - нет ошибки, 2 - нулевое расстояние м/д передатчиком и приёмником при их нулевой высоте,
			5 - нулевая длина вектора E волны, 10 - выход за пределы моделируемого интервала времени для заданного t_nTd);
		все переменные и массивы, на которые были переданы указатели в качестве входных параметров (см. описание входных параметров).
	*/


	double RG_t, RG_r, i_x, i_y, j_x, j_y, Rx, Ry, r1, r2, r0x,r0y,r0z,r,xEt,yEt,zEt,xEr,yEr,zEr,d_2,x,la,x0,y0;
	double r1x,r1y,r1z,r2x,r2y,r2z,n0x,n0y,E1,K,P_disp_, P2, rmn_x, rmn_y, r_t_mn, a1, a2, a3,Lr,k0,dd;
	double r1x0,r1y0,r1z0,nExt,nEyt,nEzt;
	unsigned long N_delta, m_, n_, N_p, N_p_max;
	unsigned long N_Np=TR.N/TR.Np; // во сколько раз количество точек в прямом и зеркально отражённом сигналах больше, чем в рассеянном сигнале (д.б. целым)
	int err_no;

	double phasa1 = 0, phasa2 = 0;
	if(n%(N_Np)==0) //условие необходимости вычисления точки рассеянного сигнала при данном n при отличии частот дискретизации для моделирования прямого и рассеянного сигналов
	{
		struct coords C;
		
		err_no=midpoint(t_nTd,&T,&R,&TR,&C);//вычисление координат приёмника и передатчика в момент t_nTd (момент излучения на n-й итерации)
		if(err_no==10) return 10;// выход сигнала за пределы моделируемого интервала
			
		double P1=P0*T.G*R.G*pow2(TR.lambda)/16./M_PI/M_PI*pow2(TR.delta)/8./M_PI; //частичные вычисления по ф-ле (2.29)
				
		RG_t=3570.*sqrt(C.zt); //радиус радиогоризонта для антенны передатчика
		RG_r=3570.*sqrt(C.zr); //радиус радиогоризонта для антенны приёмника
		
		// вычисление координат вектора r (рис. 2.9), направленного от передатчика к приёмнику:
		r0x=C.xr-C.xt;
		r0y=C.yr-C.yt;
		r0z=C.zr-C.zt;
		
		r=sqrt(pow2(r0x)+pow2(r0y)); //расстояние в проекции на плоскость Земли между передатчиком и приёмником
		if(r!=0)
		{
			Rx=r0x/r; // {Rx,Ry} - единичный вектор, направленный из точки проекции передатчика на точку проекции приёмника (проекция на плоскость Земли)
			Ry=r0y/r;
		}
		else
		{
			Rx=0.;
			Ry=1.;
		}
		
		//вычисление ортов i и j осей локальной системы координат, направленных по сторонам элементарной площадки (рис. 2.9 в 1-м отчёте):
		i_x=Ry+Rx;
		i_y=-Rx+Ry;
		r=sqrt(pow2(i_x)+pow2(i_y));
		i_x/=r;
		i_y/=r;

		j_x=-Ry+Rx;
		j_y=Rx+Ry;
		r=sqrt(pow2(j_x)+pow2(j_y));
		j_x/=r;
		j_y/=r;

		r=sqrt(pow2(r0x)+pow2(r0y)+pow2(r0z));// длина вектора r (рис. 2.9)
		if(r==0 && C.zt==0) return 2;
		if(r!=0) //если приёмник и передатчик в одной точке (r=0), то коэффициент рассогласования не считается
		{
			r0x/=r;//нормирование вектора r на 1
			r0y/=r;
			r0z/=r;
			
			// вычисление по ф-ле (2.18) координат вектора E линейно поляризованного поля антенны передатчика в системе Oxyz в направлении на приёмник:
			r1=pn(T.xp,T.yp,T.zp,r0x,r0y,r0z,&xEt,&yEt,&zEt);
			
			// вычисление по ф-ле (2.18) координат вектора E линейно поляризованного поля антенны приёмника (как если бы она была излучающей) в системе Oxyz в направлении на передатчик:
			r2=pn(R.xp,R.yp,R.zp,-r0x,-r0y,-r0z,&xEr,&yEr,&zEr);
			
			// нормировка векторов E для передатчика и приёмника, если они оба определены в направлении прямого луча (не нулевые),
			// иначе вычисление и нормировка их в направлениях зеркального луча (рис. 2.8):
			if(r1!=0 && r2!=0)
			{
				xEt/=r1;
				yEt/=r1;
				zEt/=r1;
		
				xEr/=r2;
				yEr/=r2;
				zEr/=r2;

			}
			else
			{
				d_2=pow2(C.xr-C.xt)+pow2(C.yr-C.yt); //квадрат расстояния в проекции на плоскость Земли между приёмником и передатчиком
				x=C.zt*sqrt(d_2)/(C.zt+C.zr);
				la=x/(sqrt(d_2)-x);
				x0=(C.xt+la*C.xr)/(1+la);// координаты
				y0=(C.yt+la*C.yr)/(1+la);// точки отражения зеркального луча (формулы внизу стр. 36 отчёта)
		
				r1x=x0-C.xt;// координаты вектора, направленного
				r1y=y0-C.yt;// от передатчика к точке зеркального отражения
				r1z=-C.zt;  //
				r=sqrt(pow2(r1x)+pow2(r1y)+pow2(r1z)); //длина этого вектора
				r1x/=r; //нормировка этого вектора
				r1y/=r; //
				r1z/=r; //
				
				r2x=C.xr-x0;// координаты вектора, направленного
				r2y=C.yr-y0;// от точки зеркального отражения к приёмнику
				r2z=C.zr;   //
				r=sqrt(pow2(r2x)+pow2(r2y)+pow2(r2z)); //длина этого вектора
				r2x/=r; //нормировка этого вектора
				r2y/=r; //
				r2z/=r; //

				// вычисление по ф-ле (2.18) координат вектора E линейно поляризованного поля антенны передатчика в системе Oxyz в направлении на точку зеркального отражения (x0,y0):
				r1=pn(T.xp,T.yp,T.zp,r1x,r1y,r1z,&xEt,&yEt,&zEt);
				
				if(r1==0) return 5; // нулевая длина вектора E (невозможно будет вычислить поляризацию, но это почти нереальная ситуация)
				xEt/=r1; //нормировка вектора E
				yEt/=r1; //
				zEt/=r1; //

				n0x=-r1y;// вектор нормали к вектору r1, параллельный плоскости Земли (r1={r1x,r1y} - вектор,направленный от передатчика к точке зеркального отражения)
				n0y=r1x;// (т.е. нормальный вектор плоскости падения)
				r=sqrt(pow2(n0x)+pow2(n0y));// его длина
				n0x/=r;// нормировка этого вектора
				n0y/=r;//
				E1=xEt*n0x+yEt*n0y;// горизонтальная составляющая вектора E_t={xEt,yEt} (проекция на нормальный вектор)
				xEt=2.*E1*n0x-xEt;// координаты вектора E' (рис. 2.8) (вектор поля в отражённой волне) (пишутся в переменные xEt, yEt, т.к. координаты вектора E_t больше не понадобятся)
				yEt=2.*E1*n0y-yEt;// (z-координата не меняется)
				
				// вычисление плоскости поляризации антенны приёмника (в системе Oxyz в направлении на точку зеркального отражения (x0,y0)):
				r2=pn(R.xp,R.yp,R.zp,-r2x,-r2y,-r2z,&xEr,&yEr,&zEr);
				
				if(r2==0) return 5; // нулевая длина вектора E (невозможно будет вычислить поляризацию, но это почти нереальная ситуация)
				xEr/=r2; //нормировка вектора E
				yEr/=r2; //
				zEr/=r2; //
			}

			
			// вычисление коэффициента рассогласования по амплитуде для рассеянной компоненты сигнала по ф-ле (2.22):
			K=1./sqrt(3.)+(1.-1./sqrt(3.))*fabs(xEt*xEr+yEt*yEr+zEt*zEr);
			P1*=K*K; //учёт коэффициента рассогласования
		} //if(r!=0)
				
		// количество элементарных площадок, составляющих одну сторону квадрата, вписанного
		// в окружность с радиусом радиогоризонта антенны приёмника (рис. 2.9):
		N_delta=(unsigned long)(2.*floor(RG_r/sqrt(2.)/TR.delta));// (3-я сверху формула на стр. 38)
		
		P_disp_=0.;// инициализация переменной, в которой "накапливается" суммарная мощность, рассеянная всеми площадками в направлении приёмной антенны в момент времени n+N_p_max (рассеяние сигнала, излученного передатчиком в момент времени n)
		N_p_max=0; // инициализация переменной, в которой будет максимальная задержка рассеянного луча в отсчётах

		double sin__;
		unsigned long nn;
		sin__=0.;// инициализация переменной для вычисления среднего значения суммы синусов углов падения и отражения для площадок
		nn=0;// инициализация счётчика рассеивающих площадок

		// перебор всех элементарных площадок:
		for(m_=1;m_<=N_delta;m_++)
		{
			for(n_=1;n_<=N_delta;n_++)
			{
				a1=2.*m_-1.-N_delta;
				a2=2.*n_-1.-N_delta;
				a3=TR.delta/2.;
				rmn_x = C.xr + a3*(a1*i_x+a2*j_x);// координаты центра площадки
				rmn_y = C.yr + a3*(a1*i_y+a2*j_y);// в СК Oxyz
				
				// вектор от проекции передатчика к выделенной площадке (r_mn на рис. 2.9):
				r1x=rmn_x-C.xt; a1=pow2(r1x);
				r1y=rmn_y-C.yt; a2=pow2(r1y);
				r_t_mn=sqrt(a1+a2);// длина этого вектора
								
				if(r_t_mn < RG_t) // если расстояние от передатчика до площадки < радиуса радиогоризонта передатчика (такую площадку нужно учесть):
				{
					P2=P1;// частичный расчёт формулы (2.29)
				
					r1z=-C.zt;// 3-я координата вектора r1={r1x,r1y,r1z}, направленного от передатчика к выделенной площадке
					r2x=C.xr-rmn_x;// вектор от выделенной площадки к приёмнику
					r2y=C.yr-rmn_y;//
					r2z=C.zr;	  //
					
					r1=sqrt(a1+a2+pow2(r1z));// расстояние между передатчиком и площадкой
					r2=sqrt(pow2(r2x)+pow2(r2y)+pow2(r2z));// расстояние между приёмником и площадкой
										
					Lr=r1+r2;// длина пути рассеянного площадкой луча от передатчика до приёмника
					N_p=(unsigned long)ROU((t_nTd+Lr/3e+8)/TR.Tdp); // момент прихода (в отсчётах) рассеянного площадкой 
					// луча в направлении на приёмник
					
					if(N_p<TR.Np)// если задержка луча не превышает временной интервал моделирования (не превышает размер выходного массива):
					{
						if(N_p>N_p_max) N_p_max=N_p;// максимальная задержка в отсчётах рассеянного луча (если задержка в текущей итерации цикла получилась больше, чем во всех предыдущих итерациях, то обновляем максимальное значение задержки в переменной N_p_max)
						
						P2/=pow2(r1)*pow2(r2);// частичный расчёт формулы (2.29)
						
						if(T.i==0)
						{
							dd = My_directional_diagram_3(r1x,r1y,r1z,T.teta0,T.phi0,T.dteta,T.dphi,T.n_phi,T.n_teta,T.C_phi,T.C_teta,T.C, debug);
							dd = sqrt(dd);
				//			dd=1.;
						}
						else dd=1.;// если антенна передатчика изотропная, то ДН=1
						if(dd==2) return 2; // невозможно вычислить ДН, если длина вектора r1 нулевая
						P2*=(dd*dd);// частичный расчёт формулы (2.29)
						
						if(R.i==0) {
							dd=My_directional_diagram_3(-r2x,-r2y,-r2z,R.teta0,R.phi0,R.dteta,R.dphi,R.n_phi,R.n_teta,R.C_phi,R.C_teta,R.C, debug); //вычисление ДН антенны приёмника в направлении на точку зеркального отражения, если антенна приёмника не изотропная
							My_directional_diagram_2(TR,-r2x,-r2y,-r2z, R.teta0,R.phi0, R.R_ves, R.size_ves, TR.d1, R.R_amp,
								R.R_faza, R.R_size_amp, TR.d2, &phasa2, R.MaxDN_R, chanel, mas, debug);
							dd = sqrt(dd);	
				//			dd=1.;
						}
						else dd=1.;// если антенна приёмника изотропная, то ДН=1
						if(dd==2) return 2; // невозможно вычислить ДН, если длина вектора r2 нулевая
						P2*=(dd*dd);// частичный расчёт формулы (2.29)
						*D = (dd*dd);
						sin__+=C.zt/r1+C.zr/r2; nn++;// суммирование синусов углов падения и отражения для площадок и инкремент счётчика площадок
						P2*=C.zt/r1+C.zr/r2;// частичный расчёт формулы (2.29) (sin + sin в формуле для эффективной площади рассеяния)
													
						// нормировка вектора, направленного от передатчика к выделенной площадке:
						r1x0=r1x/r1;
						r1y0=r1y/r1;
						r1z0=r1z/r1;
						
						// определение направления вектора E линейно поляризованной волны антенны передатчика в направлении на площадку:
						r=pn(T.xp,T.yp,T.zp,r1x0,r1y0,r1z0,&nExt,&nEyt,&nEzt);
													
						if(r!=0)// нормировка этого вектора, если он определён (не нулевой), в направлении на текущую площадку
						{
							nExt/=r;
							nEyt/=r;
							nEzt/=r;
							k0=1.;
						}
						else k0=0.;// если он получился нулевой в направлении на данную площадку (это м.б. в двух направлениях, т.к. приближение дипольное), то данная площадка не учитывается (k0=0, если же k0=1, то учитывается)
						
						n0x=-r1y;// вектор нормали к вектору r1, параллельный плоскости Земли
						n0y=r1x;// (т.е. нормальный вектор плоскости падения)
						r=sqrt(pow2(n0x)+pow2(n0y));// длина вектора нормали
						if(r!=0.)// если длина этого вектора не нулевая:
						{
							n0x/=r;
							n0y/=r;
							E1=nExt*n0x+nEyt*n0y;// горизонтальная составляющая вектора E_t (проекция на нормальный вектор)
							E1*=E1;// возведение в квадрат
						}
						else
						{
							E1=1.;// если длина вектора нормали нулевая, то чтобы не было неопределённости, примем горизонтальную составляющую вектора поля за 1
						}
													
						P2*=E1*TR.g0_g+(1.-E1)*TR.g0_v;// вычисление коэффициента обратного рассеяния при нормальном падении (1-я снизу ф-ла на стр. 37; на 2 здесь делить не нужно, т.к. 1/2 уже учтена в константе P1)
							
						if(s==0) P_disp[N_p]+=P2*k0;//мощность сигнала, рассеянного текущей площадкой в направлении приёмника, добавляется в ячейку массива, определяемую запаздыванием данного луча
						else P_disp_+=P2*k0;/* мощность луча, рассеянного текущей площадкой, добавляется
							 в переменную P_disp_, в которой "накапливается" суммарная мощность, рассеянная всеми
							 площадками в направлении приёмной антенны в момент времени N_p_max (рассеяние сигнала,
							 излученного передатчиком в момент времени n). Здесь делается предположение, что
							 в течение отрезка времени от n до n+N_p_max (т.е. от момента n до момента, пока луч
							 с самой большой задержкой не придёт в точку приёма) мощность сигнала, приходящего на
							 рассеивающие площадки, не меняется; поэтому в момент n+N_p_max полная рассеянная
							 мощность (от отрезка сигнала, излученного в окрестности момента n) достигает
							 наибольшего значения
							 */
						
					} // if(N_p<TR.Np)
				
				} // if(r_t_mn < RG_t)
			
			} // for(n_=1;n_<=N_delta;n_++)
		} // for(m_=1;m_<=N_delta;m_++)

		if(s==1) //если параметр s=1, то это значит, что передатчик неподвижный и рассеянная мощность принимается передатчиком
		{		 //(рассеиваемый сигнал также от передатчика), и достаточно вычислить мощность рассеянного сигнала только один раз
			for(n_=N_p_max;n_<TR.Np;n_++) // в момент N_p_max мощность установится на постоянном уровне (в этот момент придёт самый длинный луч)
			{
				P_disp[n_]=P_disp_; // установившийся уровень мощности копируется на все последующие моменты времени
				sin_mean_[n_]=sin__/nn; // вычисление среднего значения суммы синусов углов падения и отражения для площадок
			}
		}
		else sin_mean_[N_p_max]=sin__/nn;// вычисление среднего значения суммы синусов углов падения и отражения для площадок
		
		*N_delta_2=nn;

	} // if(n%(N_Np)==0)
return 0;
}
