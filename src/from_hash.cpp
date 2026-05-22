#include <cmath>
#include "CI.h"
#include <iostream>
#include <ctime>
#include <sys/time.h>
#include <iomanip>
#include <omp.h>
#include <vector>
# include "blas_link.h"
#include "molecule.h"
#include "matr.h"
#include "from_hash.h"
#include "timer.h"

# define max(a,b)  (((a)<(b))?(b):(a))
# define min(a,b)  (((a)>(b))?(b):(a))


#define HASH_SIZE 165636900UL
#define HASH_SIZE_P 15000U
typedef std::pair<int *, double *> P_T;
typedef std::pair<int, int> sign_index;

using std::cout;
using std::endl;
using std::cerr;
using std::tgammal;

double full_wall_time;
double full_cpu_time;
double wall_time;
double cpu_time;

// Считает необходимые вспомогательные факториалы. Заносит в переменную f


// Генерирует вектор размера [Na * (na + 1)], где v[i,j+1] = j заполненная орбиталь для i-ого детерминанта
// v[:, 0] = 0 удобно для реализации. Работает внутри alpha или beta подпространства
void get_vec(const int& n, const int& N, const int& size, int * v){
    if(size==1){
        v[0]=1;
//         int a=1;
//         if(n==0)a=0;
        for(int i=1;i<n+1;i++)v[i]=i;
        return;
    }
    
    
    for (long i=0; i<N*(n+1); i++) v[i]=0;
    for (int i=0; i<=n; i++) v[i]=i;
    for (long p=0; p<size-1; p++)
        for (int i=0; i<n; i++)
            if (v[ p*(n+1) + n-i] != N-i)
            {
                for (int j=1; j<n-i; j++) v[(p+1)*(n+1) + j] = v[ p*(n+1) + j];
                for (int j=n-i; j<n+1; j++) v[(p+1)*(n+1) + j] = v[ p*(n+1) + n-i] + 1 + j - (n-i);
                break;
            }
}


void hash_from_ci(const int& N, int n_empty, const int& na, const int& nb, ci_map & ci_hash, const int& i_state, double *& ci){
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1)); 
//     int * fa = new int [na * N];
//     int * fb = new int [nb * N];
//     get_factorials(na, N, fa);
//     get_factorials(nb, N, fb);
//     int * buf_a = new int [na+1];
//     int * buf_b = new int [nb+1];
//     buf_a[0]=0;
//     buf_b[0]=0;
//     ci = new double [Na*Nb];
    
//     int * fa = new int [na * N];
//     int * fb = new int [nb * N];
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    
    // Считаем факториалы
//     get_factorials(na, N, fa);
//     get_factorials(nb, N, fb);
    
    // Находим заполненные орбитали
    get_vec(na, N, Na, vec_a);
    get_vec(nb, N, Nb, vec_b);
    
    
    
    
//     for (auto i = 0u; i<Na*Nb; i++) ci[i]=0.0;
    double norm=0;
    for(int i=0;i<Na;i++)
    for(int j=0;j<Nb;j++)
    {
        ci_key key;
        key.reset();
        for(int k=1;k<=na;k++)
            key[vec_a[i*(na+1)+k]-1+n_empty]=1;
        
        for(int l=1;l<=nb;l++)
            key[vec_b[j*(nb+1)+l]-1+n_empty+CI_MAX_SPACE]=1;
            
//         ci_hash[key] = ci[i*Nb+j];
        ci_hash[key] = ci[i*Nb+j];
//         norm+=ci[i*Nb+j]*ci[i*Nb+j];
    }
//     printf("S_arr = %.10e\n",norm);
//     norm=0;
    for (const auto& [key, value] : ci_hash){
        norm+=value*value;
    }
    
//     printf("S_hash = %.10e\n",norm);

    
    
    
    
//    delete[] fa    ;
//    delete[] fb    ;
   delete[] vec_a ;
   delete[] vec_b ;
}

void occ_from_ci(const int& N, const int& na, const int& nb, std::vector<std::vector<int >> * occ){
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1)); 
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
 
    get_vec(na, N, Na, vec_a);
    get_vec(nb, N, Nb, vec_b);
    
    double norm=0;
    for(int i=0;i<Na;i++)
    for(int j=0;j<Nb;j++)
    {
        for(int k=0;k<2*N;k++)occ[0][i*Nb+j][k]=0;
        for(int k=1;k<=na;k++)
            occ[0][i*Nb+j][vec_a[i*(na+1)+k]-1]=1;
        
        for(int l=1;l<=nb;l++)
            occ[0][i*Nb+j][vec_b[j*(nb+1)+l]-1+N]=1;
        
    }

   delete[] vec_a ;
   delete[] vec_b ;
}

void occ_from_ci_ba(const int& N, const int& na, const int& nb, std::vector<std::vector<int >> * occ){
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1)); 
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    
    //fprintf(stderr, "start index generation for na = %d, nb = %d (inverted)\n",na,nb);
    
    get_vec(na, N, Na, vec_a);
    get_vec(nb, N, Nb, vec_b);
    
    double norm=0;
    for(int i=0;i<Nb;i++)
    for(int j=0;j<Na;j++)
    {
        for(int k=0;k<2*N;k++)occ[0][i*Na+j][k]=0;
        for(int k=1;k<=na;k++)
            occ[0][i*Na+j][vec_a[j*(na+1)+k]-1]=1;
        
        for(int l=1;l<=nb;l++)
            occ[0][i*Na+j][vec_b[i*(nb+1)+l]-1+N]=1;
    }

   delete[] vec_a ;
   delete[] vec_b ;
}


void ci_from_hash(const int& N, const int& na, const int& nb, ci_map_arr& ci_hash, const int& i_state, double *& ci){
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1)); 
    int * fa = new int [na * N];
    int * fb = new int [nb * N];
    get_factorials(na, N, fa);
    get_factorials(nb, N, fb);
    int * buf_a = new int [na+1];
    int * buf_b = new int [nb+1];
    buf_a[0]=0;
    buf_b[0]=0;
    ci = new double [Na*Nb];
    for (auto i = 0u; i<Na*Nb; i++) ci[i]=0.0;
    for (const auto& [key, value] : ci_hash){
        auto ind_a = get_index(key, na, N, fa, buf_a, 0);
        auto ind_b = get_index(key, nb, N, fb, buf_b, CI_MAX_SPACE);
        ci[ind_a*Nb + ind_b] = value[i_state];
    } 
    
    
    
    delete[] fa;
    delete[] fb;
    delete[] buf_a;
    delete[] buf_b;
}

// Функции для подсчета времени. Wall time и cpu time
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)std::clock() / CLOCKS_PER_SEC;
}

// Start time инициализирует счетчики
void start_time(){
    full_wall_time = get_wall_time();
    full_cpu_time = get_cpu_time();
    wall_time = get_wall_time();;
    cpu_time = get_cpu_time();
}

// Вывод в стрим времени прошедшего с вызова start_time (TOTAL)
// и с предыдущего вызова time (ON STEP)
void time(std::ostream& os, const char * description){
    double wall_buf = get_wall_time();
    double cpu_buf = get_cpu_time();
    os<<"--------"<<description<<"--------"<<endl;
    os<<"  TOTAL ->  WALL:"<<std::setw(10)<<std::fixed<<std::setprecision(2)<<wall_buf - full_wall_time<<"s.    CPU:"<<std::setw(10)<<std::fixed<<std::setprecision(2)<<cpu_buf - full_cpu_time<<"s."<<std::endl;
    os<<"ON STEP ->  WALL:"<<std::setw(10)<<std::fixed<<std::setprecision(2)<<wall_buf - wall_time<<"s.    CPU:"<<std::setw(10)<<std::fixed<<std::setprecision(2)<<cpu_buf - cpu_time<<"s."<<std::endl;
    wall_time = wall_buf;
    cpu_time = cpu_buf;
}



// Расчет Г коэффициента Г_i = П_{q<i} (-1)**(k_q)
// Минус один в степени количества занятых орбиталей меньше k


// Считает вклад в гамильтониан от диагональных, однократно- и двукратно- возбужденных alpha-alpha (или beta-beta) орбиталей.
// Буквой s обозначается интересующее нас спиновое подпространство.
// Буквой p - второе спиновое подпространство. 
// В частности ci_pp - скалярное произведение КВ коэффициентов по всем детерминантам p пространства для некоторых двух детерминантов s пространства.
// Диагональный и однократно-возбужденные вклады считаются только внутри орбиталей заданного спинового подпространства. 
// ВСЕ alpha-beta вклады выносятся в функцию H_sp
// inc - шаг между КВ коэффициентами p подпространства внутри матрицы хранения для BLAS. inc_alpha = 1, inc_beta = Na.
// mult - множитель сдвига между КВ коэффициентами s подпространства внутри матрицы хранения для BLAS. mult_alpha = Na, mult_beta = 1.
// Ns и Np - размеры подпространств s и p.
double H_ss(double * H, const int& N, const int& ns, const int& Ns, int * vec_s, int * fs, const int& Np,
               const double * ci_1, const double * ci_2, const double * G, int n_s1, int ld1, int n_s2, int ld2, int nt, double c){
    // Инициализируем вспомогательные переменные
    int N3 = N*N*N;
    int N2 = N*N;
    double result = 0.0;
    double ** r;
    r = new double *[nt];
    for (int i=0; i<nt; i++) r[i]=new double[n_s1*n_s2];
    for (int i=0; i<nt; i++) set_zero_matr(r[i],n_s1*n_s2);
    // Проходимся по всем s векторам
    #pragma omp parallel for
    for (int i = 0; i<Ns; i++){
//         printf("th = %d\n",omp_get_thread_num());
        
        double * H_th = r[omp_get_thread_num()];
        double * ci_pp;
        ci_pp = new double[n_s1*n_s2];
        int * buf = new int [ns + 1];
        int * vec = new int [N];
        buf[0]=0;
        // Генерируем ON вектор по заселенным орбиталям i детерминанта.
//         for (int vec_ind=0; vec_ind<N; vec_ind++) vec[vec_ind] = 0;
        memset(vec,0,N*sizeof(int));
        for (int k=1; k<ns+1; k++) vec[vec_s[i*(ns+1)+k]-1] = 1;
        
        // Вклад ss орбиталей в диагональный элемент.
            // Сумма pp ci коэффициентов для диагонального ii детерминанта
//             double ci_pp = cblas_ddot( Np, ci_1+i*Np, 1 , ci_2+i*Np , 1);
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        n_s1,n_s2,Np,1.0,
                        ci_1+i*Np*ld1,ld1,
                        ci_2+i*Np*ld2,ld2,0.0,
                        ci_pp,n_s2);
            
            // Два цикла по орбиталям ON вектора i
            for (int k=0; k<N; k++)
                for (int l=0; l<N; l++){
                    // Не используем if конструкции. Домножаем на заселенности.
                    double tmp = vec[k]*vec[l]*0.5*(G[k*N3+k*N2+l*N+l] - G[k*N3+l*N2+l*N+k]);
                    for(int i_s=0;i_s<n_s1;i_s++)
                    for(int j_s=0;j_s<n_s2;j_s++)
                    H_th[i_s*n_s2+j_s] += tmp*ci_pp[i_s*n_s2+j_s];
                }
        // Однократные и двукратные возбуждения
            for (int k=0; k<N; k++)
                if (vec[k]){
                    // Получаем Г_i^k множитель
                    auto G_k = get_sign_from_ON(vec, k);
                    // Проходимся по орбиталям и ищем пустые для однократного возбуждения
                    for (int l=0; l<N; l++)
                        if (!(vec[l])){
                            // Находим индекс орбитали и Г_ind_l^l множитель
                            vec[k] = 0;
                            vec[l] = 1;
                            auto G_l = get_sign_from_ON(vec, l);
                            auto ind_l = get_ind_from_ON(vec, N, ns, fs,  buf);
                            // Сумма pp ci коэффициентов для недиагонального i ind_l детерминанта
//                             ci_pp = cblas_ddot( Np, ci_1+i*Np , 1 , ci_2+ind_l*Np , 1);
                            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                                        n_s1,n_s2,Np,1.0,
                                        ci_1+i    *Np*ld1,ld1,
                                        ci_2+ind_l*Np*ld2,ld2,0.0,
                                        ci_pp,n_s2);
            
                            // Вклад от однократных возбуждений по ss орбиталям. Еще раз цикл по орбиталям. 
                            // Не используем if конструкции. Домножаем на заселенность j орбитали.
                            // Флипаем l орбиталь. l и k в vec теперь пустые
                            vec[l] = 0;
                            for (int j=0; j<N; j++){
                                double tmp = vec[j]*G_k*G_l*(G[k*N3+l*N2+j*N+j] - G[k*N3+j*N2+j*N+l]);
                                for(int i_s=0;i_s<n_s1;i_s++)
                                for(int j_s=0;j_s<n_s2;j_s++)
                                    H_th[i_s*n_s2+j_s] += tmp*ci_pp[i_s*n_s2+j_s];
                            }
                            // Начинаем двукратные возбуждения
                            // Цикл только по орбиталям большим, чем k
                            for (int kk=k+1; kk<N; kk++)
                                if (vec[kk]){
                                    // Получаем Г_i^kk множитель.
                                    vec[k] = 1; //1
                                    auto G_kk = get_sign_from_ON(vec, kk);
                                    vec[l] = 1; //1
                                    // Цикл только по большим, чем l
                                    for (int ll=l+1; ll<N; ll++)
                                        if (!(vec[ll])){
                                            // Находим индекс орбитали и Г_ind_ll^ll множитель. Г_ind_ll^l != Г_ind_l^l, его  нужно считать
                                            vec[k ] = 0;
                                            vec[kk] = 0;
                                            vec[ll] = 1;
                                            // k=0, kk=0, ll=1, l=1
                                            G_l = get_sign_from_ON(vec, l);
                                            auto G_ll = get_sign_from_ON(vec, ll);
                                            auto ind_ll = get_ind_from_ON(vec, N, ns, fs,  buf);
                                            // Сумма pp ci коэффициентов для недиагонального i ind_ll детерминанта
//                                             ci_pp = cblas_ddot( Np, ci_1+i*Np , 1 , ci_2+ind_ll*Np , 1);
                                            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                                                        n_s1,n_s2,Np,1.0,
                                                        ci_1+i     *Np*ld1,ld1,
                                                        ci_2+ind_ll*Np*ld2,ld2,0.0,
                                                        ci_pp,n_s2);
            
                                            // Вклад от двукратного возбуждения.
//                                             result += G_k*G_kk*G_l*G_ll*(G[k*N3+l*N2+kk*N+ll] - G[kk*N3+l*N2+k*N+ll])*ci_pp;
                                            double tmp = G_k*G_kk*G_l*G_ll*(G[k*N3+l*N2+kk*N+ll] - G[kk*N3+l*N2+k*N+ll]);
                                            for(int i_s=0;i_s<n_s1;i_s++)
                                            for(int j_s=0;j_s<n_s2;j_s++)
                                                H_th[i_s*n_s2+j_s] += tmp*ci_pp[i_s*n_s2+j_s];
                                            // Возвращаем вектор с l=1
                                            vec[k ] = 1;
                                            vec[kk] = 1;
                                            vec[ll] = 0;
                                            // k=1, kk=1, ll=0, l=1
                                        }
                                    // Возвращаем l=0;
                                    vec[l] = 0; //0
                                }
                            vec[k] = 1;
                        }
                }
        delete[] buf;
        delete[] vec;
        delete[] ci_pp;
    }
    for(int i=0; i<nt; i++)
    for(int i_s=0;i_s<n_s1;i_s++)
    for(int j_s=0;j_s<n_s2;j_s++)
        H[i_s*n_s2+j_s]+=r[i][i_s*n_s2+j_s]*c;
    
    for (int i=0; i<nt; i++) delete[] r[i];
    delete[] r;

    
    return result;
}


double H_sp(double * H,
            const int& N, 
            const int& na,  const int& Na, int * vec_a, int * fa, 
            const int& nb,  const int& Nb, int * vec_b, int * fb, 
            double * ci_1, double * ci_2, double * G, int n_s1, int ld1, int n_s2, int ld2, int nt){
    // Вспомогательные переменные.
    const int Na_zv = (int) std::lround(tgammal(N) / tgammal(na) / tgammal(N-na+1));
    const int Nb_zv = (int) std::lround(tgammal(N) / tgammal(nb) / tgammal(N-nb+1));
    const int N3 = N*N*N;
    const int N2 = N*N;
    int * Ia = new int [Na_zv];
    int * Ia_friends = new int [Na_zv * (N-na) * 3];
    int * Ib = new int [Nb_zv];
    int * Ib_friends = new int [Nb_zv * (N-nb) * 3];
    double result = 0.0;
    double ** r;
    r = new double *[nt];
    for (int i=0; i<nt; i++) r[i]=new double[n_s1*n_s2];
    for (int i=0; i<nt; i++) set_zero_matr(r[i],n_s1*n_s2);
    
    //cerr<<endl<<"--------Start H_sp--------"<<endl;
    // Цикл по alpha орбиталям
    for (int i=0u; i<N; i++){
        // Для i-ой alpha орбитали создаем массив занятых alpha детерминантов и массив однократных возбуждений
        get_I_I_friends(na, N, Na, fa, i, vec_a, Ia, Ia_friends);
        // Цикл по beta орбиталям
        for (int j=0u; j<N; j++){
            // Для j-ой beta орбитали создаем массив занятых beta детерминантов и массив однократных возбуждений
            get_I_I_friends(nb, N, Nb, fb, j, vec_b, Ib, Ib_friends);
            // Проходимся по массивам alpha и beta детерминантов с занятыми i и j орбиталями соотвестсвенно
            for (int k = 0; k<Na_zv; k++){
                double result_omp = 0.0;
                #pragma omp parallel for 
                for (int l = 0; l<Nb_zv; l++){
                    double * H_th = r[omp_get_thread_num()];

                    // Переменные соотвующие индексам КВ массива
//                     printf("k = %d, l = %d\n",k,l);
                    int ind_k = Ia[k], ind_l = Ib[l];
//                     printf("k = %d, l = %d, ind_k = %d, Ia[k] = %d, ind_l = %d, Ib[l] = %d\n", k ,l, ind_k, Ia[k], ind_l, Ib[l]);
//                     getchar();
                    // Вклад в диагональный элемент по i и j орбиталям от ind_k ind_l детерминанта
                    for(int i_s=0;i_s<n_s1;i_s++)
                    for(int j_s=0;j_s<n_s2;j_s++)
                        H_th[i_s*n_s2+j_s] += ci_1[(ind_k*Nb + ind_l)*ld1+i_s ] * 
                                              ci_2[(ind_k*Nb + ind_l)*ld2+j_s ] * G[ i*N3 + i*N2 + j*N + j ];
                    
                    // Вклад в H по однократным возбуждениям j
                    for (int n = 0; n<(N-nb); n++){
                        // ind_n указывает на 3 числа: индекс КВ возбужденного отн l детерминанта по j beta орбитали, орбиталь возбуждения, знак перехода.
                        int * ind_n{Ib_friends + l*(N-nb)*3 + n*3};
                        for(int i_s=0;i_s<n_s1;i_s++)
                        for(int j_s=0;j_s<n_s2;j_s++)
                            H_th[i_s*n_s2+j_s] += ci_1[(ind_k*Nb + ind_l   )*ld1+i_s]*
                                                  ci_2[(ind_k*Nb + ind_n[0])*ld2+j_s]*G[i*N3+i*N2+j*N+ind_n[1]]*ind_n[2]; 
                    }
                    
                    // Вклад в H по однократным возбуждениям i и по двукратным возбуждениям i и j
                    for (int m = 0; m<(N-na); m++){
                        // ind_m указывает на 3 числа: индекс КВ возбужденного по k детерминанта по i орбитали, орбиталь возбуждения, знак перехода.
                        int * ind_m{Ia_friends + k*(N-na)*3 + m*3};
                        for(int i_s=0;i_s<n_s1;i_s++)
                        for(int j_s=0;j_s<n_s2;j_s++)
                            H_th[i_s*n_s2+j_s] += ci_1[(ind_k   *Nb + ind_l)*ld1+i_s]*
                                                  ci_2[(ind_m[0]*Nb + ind_l)*ld2+j_s]*G[j*N3+j*N2+i*N+ind_m[1]]*ind_m[2];
                        for (int n = 0; n<(N-nb); n++){
                            // ind_n указывает на 3 числа: индекс КВ возбужденного по l детерминанта по j орбитали, орбиталь возбуждения, знак перехода.
                            int * ind_n{Ib_friends + l*(N-nb)*3 + n*3};
                            for(int i_s=0;i_s<n_s1;i_s++)
                            for(int j_s=0;j_s<n_s2;j_s++)
                                H_th[i_s*n_s2+j_s] += ci_1[(ind_k   *Nb + ind_l   )*ld1+i_s]*
                                                      ci_2[(ind_m[0]*Nb + ind_n[0])*ld2+j_s]*
                                                      G[i*N3+ind_m[1]*N2+j*N+ind_n[1]]*ind_m[2]*ind_n[2];
                        }
                    }
                }
//                 result+=result_omp;
            }
            //cerr<<"     finished "<<std::setw(6)<<std::setprecision(2)<<std::fixed<< 100.0*(i*N+j)/(N-1)/(N+1)<<"%        \r\b";
        }
    }
    for (int i=0; i<nt; i++)
    for(int i_s=0;i_s<n_s1;i_s++)
    for(int j_s=0;j_s<n_s2;j_s++)
        H[i_s*n_s2+j_s]+=r[i][i_s*n_s2+j_s];

    for (int i=0; i<nt; i++)delete[] r[i];
    delete[] r;
    
    //cerr<<endl;
    // Чистим память
    delete[] Ia;
    delete[] Ia_friends;
    delete[] Ib;
    delete[] Ib_friends; 
    
    // Возращаем result
    return result;
}


 void transpose(double * ci_1, double * ci_2, const int& Na, const int& Nb){
     if (Na==Nb){
        #pragma omp parallel for
        for (auto i = 0u; i<Na-1; i++){
            double buf = 0.0;
            for (auto j = i+1; j<Na; j++){
                buf = ci_1[i*Na+j];
                ci_1[i*Na+j] = ci_1[j*Na+i];
                ci_1[j*Na+i] = buf;
                buf = ci_2[i*Na+j];
                ci_2[i*Na+j] = ci_2[j*Na+i];
                ci_2[j*Na+i] = buf;
            }
        }
     }
     else {
        int size = Na*Nb - 1;
        double temp_1;
        double temp_2;
        int next;
        int cycleBegin;
        int i;
        std::vector<std::bitset<HASH_SIZE>> wrapper(1);
        auto & b = wrapper[0];
        b.reset(); 
        b[0] = b[size] = 1; 
        i = 1;
        while (i < size) 
        { 
            cycleBegin = i; 
            temp_1 = ci_1[i];
            temp_2 = ci_2[i];
            do
            { 
                next = (((long)i*Na)%size);
                std::swap(ci_1[next], temp_1); 
                std::swap(ci_2[next], temp_2);
                b[i] = 1; 
                i = next; 
            } 
            while (i != cycleBegin); 
            for (i = 1; i < size && b[i]; i++) ; 
        } 
    }

}


// Транспонирование одной матрицы. Надо переписать
void transpose_one_matrix(double * A, const int& Na, const int& Nb){
     if (Na==Nb){
        #pragma omp parallel for
        for (auto i = 0u; i<Na-1; i++){
            double buf = 0.0;
            for (auto j = i+1; j<Na; j++){
                buf = A[i*Na+j];
                A[i*Na+j] = A[j*Na+i];
                A[j*Na+i] = buf;
            }
        }
     }
     else {
        int size = Na*Nb - 1;
        double temp;
        int next;
        int cycleBegin;
        int i;
        std::vector<std::bitset<HASH_SIZE>> wrapper(1);
        auto & b = wrapper[0];
        b.reset(); 
        b[0] = b[size] = 1; 
        i = 1;
        while (i < size) 
        { 
            cycleBegin = i; 
            temp = A[i];
            do
            { 
                next = (((long)i*Na)%size);
                std::swap(A[next], temp); 
                b[i] = 1; 
                i = next; 
            } 
            while (i != cycleBegin); 
            for (i = 1; i < size && b[i]; i++) ; 
        } 
    }

}

double calc_2e_ints(double * H, const int& N, const int& na, const int& nb, double * ci_1, double * ci_2, double * G,
                    const bool& singlet, const int n_s1, const int ld1, const int n_s2, const int ld2){
    // Старт таймера
    start_time();
    
    // Вспомогательные переменные
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1));    
    int * fa = new int [na * N];
    int * fb = new int [nb * N];
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    
    // Считаем факториалы
    get_factorials(na, N, fa);
    get_factorials(nb, N, fb);
    
    // Находим заполненные орбитали
    get_vec(na, N, Na, vec_a);
    get_vec(nb, N, Nb, vec_b);
        
    
    double H_aa{0.0}, H_bb{0.0}, H_ab{0.0};
    
    // Расчет alpha-alpha вклада
    set_zero_matr(H,n_s1*n_s2);

    H_ss(H, N, na, Na, vec_a, fa, Nb, ci_1, ci_2, G, n_s1, ld1, n_s2, ld2, num_threads,1.0);
    // Таймер
//     time(std::cerr, "H_aa calc");
//     cerr<<"V_aa = "<<std::scientific<<std::setw(15)<<std::setprecision(10)<<H_aa<<endl<<endl;
//     exit(0);
    double * ci_1_tr = new double[Na*Nb*ld1];
    double * ci_2_tr = new double[Na*Nb*ld2];
    transpose_3d_abc_to_bac(ci_1_tr, ci_1, Na, Nb, ld1);
    transpose_3d_abc_to_bac(ci_2_tr, ci_2, Na, Nb, ld2);
    
    H_ss(H, N, nb, Nb, vec_b, fb, Na, ci_1_tr, ci_2_tr, G, n_s1, ld1, n_s2, ld2, num_threads,1.0);
    // Таймер
//     time(std::cerr, "H_bb calc");
    
    
    // Расчет alpha-beta вклада
    H_ab = H_sp(H, N, na, Na, vec_a, fa, nb, Nb, vec_b, fb, ci_1, ci_2, G, n_s1, ld1, n_s2, ld2, num_threads);
    // Таймер
//     time(std::cerr, "H_ab calc");    
//     
    // Чистим память
    delete[] fa;
    delete[] fb;
    delete[] vec_a;
    delete[] vec_b;
    delete[] ci_1_tr;
    delete[] ci_2_tr;
    
    
    return 0;
}

int gen_1_el_DMA_aldet(double * O,double * ci_A, double * ci_B, int nA, int nB, int n_act, int na, int nb){
    
    int Na = (int) std::lround(tgammal(n_act+1) / tgammal(na+1) / tgammal(n_act-na+1));
    int Nb = (int) std::lround(tgammal(n_act+1) / tgammal(nb+1) / tgammal(n_act-nb+1));    
    int * fa = new int [na * n_act];
    int * fb = new int [nb * n_act];
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    get_factorials(na, n_act, fa);
    get_factorials(nb, n_act, fb);
    
    get_vec(na, n_act, Na, vec_a);
    get_vec(nb, n_act, Nb, vec_b);
    
    int * bit_a = new int [n_act];
    int * bit_b = new int [n_act];
    int * buf = new int [max(na+1,nb+1) + 1];//more space for +a-b and -a+b ????
    buf[0]=0;
    

    double * coef_A;
    double * coef_B;
    double sign_A;
    double sign_B;
    
    double * BUF = new double[nA*nB*n_act*n_act];
    set_zero_matr(BUF,nA*nB*n_act*n_act);
    double * BUF1;
    
    
    for(int i_CI = 0; i_CI<Na; i_CI++){
        memset(bit_a,0,n_act*sizeof(int));
        for (int k=1; k<na+1; k++) bit_a[vec_a[i_CI*(na+1)+k]-1] = 1;
        
        
        sign_A=1;
        coef_A=ci_A+i_CI*Nb*nA;
        for(int i=0;i<n_act;i++)
            if(bit_a[i]){
                bit_a[i]=0;
                sign_B=1;
                    for(int j=0;j<n_act;j++)
                    if(bit_a[j]==0){
                        BUF1 = BUF + (i*n_act+j)*nA*nB;
                        bit_a[j]=1;
                        
                        auto i_CIext = get_ind_from_ON(bit_a, n_act, na, fa, buf);
                        coef_B=ci_B+i_CIext*Nb*nB;
                        
//                         for(int i_s=0;i_s<nA;i_s++)
//                         for(int j_s=0;j_s<nB;j_s++)
//                         for(int j_CI = 0; j_CI<Nb; j_CI++)
//                             BUF1[i_s*nB+j_s]+=coef_A[j_CI*nA+i_s]*coef_B[j_CI*nB+j_s]*sign_A*sign_B;
                        
                        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                                    nA,nB,Nb,sign_A*sign_B,
                                    coef_A,nA,
                                    coef_B,nB,1.0,
                                    BUF1,nB);
                        
                        
                        bit_a[j]=0;
                    }
                    else
                        sign_B=-sign_B;

                sign_A=-sign_A;
                bit_a[i]=1;
            }
    }
//     printf_timer("transpose");
    
    for(int i_s=0;i_s<nA*nB;i_s++)
    for(int i=0;i<n_act*n_act;i++)
//     for(int j=0;j<n_act;j++)
//     for(int j_s=0;j_s<nB;j_s++)
        O[i_s*n_act*n_act+i]=BUF[i*nA*nB+i_s];
                        
    
    
   
    return 0;
}

int gen_1_el_DMB_aldet(double * O,double * ci_A, double * ci_B, int nA, int nB, int n_act, int na, int nb){
    
    int Na = (int) std::lround(tgammal(n_act+1) / tgammal(na+1) / tgammal(n_act-na+1));
    int Nb = (int) std::lround(tgammal(n_act+1) / tgammal(nb+1) / tgammal(n_act-nb+1));    
    int * fa = new int [na * n_act];
    int * fb = new int [nb * n_act];
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    get_factorials(na, n_act, fa);
    get_factorials(nb, n_act, fb);
    
    get_vec(na, n_act, Na, vec_a);
    get_vec(nb, n_act, Nb, vec_b);
    
    int * bit_a = new int [n_act];
    int * bit_b = new int [n_act];
    int * buf = new int [max(na+1,nb+1) + 1];//more space for +a-b and -a+b ????
    buf[0]=0;
    

    double * coef_A;
    double * coef_B;
    double sign_A;
    double sign_B;
    
    set_zero_matr(O,nA*nB*n_act*n_act);
    
    for(int j_CI = 0; j_CI<Nb; j_CI++){
        memset(bit_a,0,n_act*sizeof(int));
        memset(bit_b,0,n_act*sizeof(int));
        for (int k=1; k<nb+1; k++) bit_b[vec_b[j_CI*(nb+1)+k]-1] = 1;
        
        sign_A=1;
        coef_A=ci_A+j_CI*nA;
        for(int i=0;i<n_act;i++)
            if(bit_b[i]){
                bit_b[i]=0;
                sign_B=1;
                    for(int j=0;j<n_act;j++)
                    if(bit_b[j]==0){
                        bit_b[j]=1;
                        
                        auto j_CIext = get_ind_from_ON(bit_b, n_act, na, fa, buf);
                        coef_B=ci_B+j_CIext*nB;
                        
                        for(int i_s=0;i_s<nA;i_s++)
                        for(int j_s=0;j_s<nB;j_s++)
                        for(int i_CI = 0; i_CI<Na; i_CI++)
                            O[((i_s*nB+j_s)*n_act+i)*n_act+j]+=coef_A[i_s+i_CI*Nb*nA]*coef_B[j_s+i_CI*Nb*nB]*sign_A*sign_B;
                        
                        bit_b[j]=0;
                    }
                    else
                        sign_B=-sign_B;

                sign_A=-sign_A;
                bit_b[i]=1;
            }
    }
   
    return 0;
}




// Расчет количества инверсий в перестановке вектора заселенных орбиталей
int inversion_sign(const int * new_vec, const int& n){
    int count = 0;
    for (int i = 1; i<=n; i++)
        for (int j = i+1; j<=n; j++)
            if (new_vec[i]>new_vec[j]) count++;
    return (count%2 == 0 ? 1 : -1);
}

// typedef std::pair<int, int> sign_index;
// Получаем знак и индекс для перестаервки орбиталей по матрице P
sign_index get_sign_index(const int& N, const int& n, int * f,  int * vec, int * new_vec, int * ON, int * P){
    new_vec[0]=0;
    for (int i=1; i<=n; i++) new_vec[i] = P[vec[i]-1]+1;
    for (int i=0; i<N; i++) ON[i]=0;
    for (int i=1; i<=n; i++) ON[new_vec[i]-1]=1;
    int sign = inversion_sign(new_vec, n);
    int index = get_ind_from_ON(ON, N, n, f, new_vec);
    return {sign, index};
}

// LU Lapack дает пары инверсий
// Нам для работы нужно немного по-другому
// Пересчитывем свое P с блэкджеком
// В данном виде P[i] указывает на столбец в матриице перестановок
// с ненулевым элементом в строке i.
void get_nice_P(const int& N, int * P){
    int * P_buf = new int[N];
    for (int i = 0; i<N; i++) P_buf[i]=i;
    for (int i = 0; i<N; i++) std::swap(P_buf[P[i]-1],P_buf[i]);
    for (int i = 0; i<N; i++) P[P_buf[i]]=i;
    delete[] P_buf;
}


// Обертка над LU Lapack
void LU(const int& size, double * A, int * P){
    int M = size;
    int N = size;
    int LDA = size;
    int INFO;
    // Переводим в column major
    transpose_one_matrix(A, size, size);
    // Lapack dirty work
    dgetf2_(&M, &N, A, &LDA, P, &INFO);
    // Транспонируем LU матрицу
    transpose_one_matrix(A, size, size);
    // Получаем нужную нам для работу
    // "правильную" P матрицу
    get_nice_P(size, P);
    for (int i =0; i<size; i++)
//         std::cerr<<P[i]<<" ";
//     std::cerr<<std::endl;
    if (INFO!=0){
        std::cerr<<"Error on LU. Ooops."<<std::endl;
        throw "Error on LU. Ooops.";
    }
    if (INFO==0){
//         std::cerr<<"Ok with LU"<<std::endl;
    }
}

// Обращение верхнетреугольной матрицы.
void InvU(const int& size, double * U){
    int * P_eq = new int[size];
    for (int i = 1; i<=size; i++) P_eq[i-1]=i;
    int N = size;
    int LDA = size;
    int LWORK = 16*size;
    int INFO;
    double * WORK = new double[LWORK];
    // Переходим в fortran style
    transpose_one_matrix(U, size, size);
    dgetri_( &N, U, &LDA, P_eq, WORK, &LWORK, &INFO);
    // Перевод результата в C-style
    transpose_one_matrix(U, size, size);
    if (INFO!=0){
        std::cerr<<"Error on InvU. Ooops."<<std::endl;
        throw "Error on InvU. Ooops.";
    }
    delete[] P_eq;
    delete[] WORK;
}



// typedef std::pair< int *, double *> P_T;
// Из матрицы получаем перехода в новые МО получаем 
// P и T матрицы для пересчета КВ коэффициентов
// Не забываем чистить память вне функции
P_T calc_P_T(const int& size, double * A){
    int * P = new int[size];
    double * T = new double[size*size];
    LU(size, A, P);
//     printf("A\n");
//     PrintMatr(A,size,size,1);
    for (int i = 0; i<size; i++)
        for (int j = 0; j<i; j++){
            T[i*size+j]=-1.0*A[i*size+j];
            A[i*size+j]=0.0;
        }
    
    InvU(size, A);
    for (int i = 0; i<size; i++)
        for (int j = i; j<size; j++) T[i*size+j]=A[i*size+j];
//     printf("T\n");
//     PrintMatr(T,size,size,1);
    return {P, T};
}


// Переставлять КВ по матрице P будем по следующему принципу.
// Проходимся по alpha, транспонируем, проходимся по beta
// Внутри прохода:
//     Создаем bitset b переставленных строк (изначально заполнен нулями)
//     Проходимся по всем строкам i. Если переставлена идем дальше (i+1)
//     Иначе начинаем новую цепочку: находим новое положение i_new и знак перехода. b(i_new)=1; 
//     Заносим строку i_new в буфер. На ее место ставим строку i.
//     Продолжаем цепочку по i_new. По окончании цепочки возвращаемся к i+1

void change_ci_with_P_one_dim(double * ci, const int& N, const int& ns, const int& Ns, int * vec_s, int * fs, const int& Np, int * P){
    double * buf_row = new double[Np];
    int * new_vec = new int[ns+1];
    int * ON = new int[N];
    int last, next, sign;
    sign_index si;
    std::vector<std::bitset<HASH_SIZE_P>> wrapper(1);
    auto & b = wrapper[0];
    b.reset();
    for (int i=0; i<Ns; i++){
        if (!b[i]){
            si = get_sign_index(N, ns, fs, vec_s+i*(ns+1), new_vec, ON, P);
            sign = si.first;
            next = si.second;
	    if ((next == i) && (sign==-1)){
	    	cblas_dscal(Np, -1.0, ci+i*Np, 1);
		continue;
	    }
            if (next != i){
                last = next;
                cblas_dswap (Np, buf_row, 1, ci+i*Np, 1);
                do{
                    if (sign==-1) cblas_dscal (Np, -1.0, buf_row, 1);
                    cblas_dswap (Np, ci+next*Np, 1, buf_row, 1);
                    b[next]=1;
                    si = get_sign_index(N, ns, fs, vec_s+next*(ns+1), new_vec, ON, P);
                    sign = si.first;
                    next = si.second;
                } while (next!=last);
            }
        }
    }
    delete[] buf_row;
    delete[] new_vec;
    delete[] ON;
}

// Функция, объеденяющая проход P по alpha и beta.
// !!! В конце ci_1 и ci_2 транспонированы
void change_ci_with_P(double * ci_1, double * ci_2, const int& N, const int& na, const int& Na, int * vec_a, int * fa, const int& nb, const int& Nb, int * vec_b, int * fb, int * P_1, int * P_2){
    change_ci_with_P_one_dim(ci_1, N, na, Na, vec_a, fa, Nb, P_1);
    change_ci_with_P_one_dim(ci_2, N, na, Na, vec_a, fa, Nb, P_2);
    transpose(ci_1, ci_2, Na, Nb);
    change_ci_with_P_one_dim(ci_1, N, nb, Nb, vec_b, fb, Na, P_1);
    change_ci_with_P_one_dim(ci_2, N, nb, Nb, vec_b, fb, Na, P_2);
}

void change_ci_with_P_1(double * ci_1, const int& N, const int& na, const int& Na, int * vec_a, int * fa, const int& nb, const int& Nb, int * vec_b, int * fb, int * P_1){
    change_ci_with_P_one_dim(ci_1, N, na, Na, vec_a, fa, Nb, P_1);
    transpose_one_matrix(ci_1, Na, Nb);
    change_ci_with_P_one_dim(ci_1, N, nb, Nb, vec_b, fb, Na, P_1);
}

// Переход к новым КВ коэффициентам по матрице T вдоль одного спинового пространства.
void change_ci_with_T_one_dim(const int& N, const int& ns, const int& Ns, int * vec_s, int * fs, const int& Np, double *& ci, double *& ci_new, double * T){
    int * buf = new int [(ns+1)];
    buf[0]=0;
    int * vec = new int [N];
    get_vec(ns, N, Ns, vec_s);
    get_factorials(ns, N, fs);
    for (int k = 0; k<N; k++){
        memset(ci_new, 0, sizeof(double)*Ns*Np);
        for (int i = 0; i<Ns; i++){
            memset(vec, 0, sizeof(int)*N);
            for (int l=1; l<ns+1; l++) vec[vec_s[i*(ns+1)+l]-1] = 1;
            double * ci_i = ci+i*Np;
            if (!vec[k]) {
                cblas_daxpy(Np, 1.0, ci_i, 1, ci_new+i*Np, 1);
                continue;
            }
            auto G_k = get_sign_from_ON(vec, k);
            cblas_daxpy(Np, T[k*N+k], ci_i, 1, ci_new+i*Np, 1);
            for (int l = 0; l<N; l++){
                if (!vec[l]){
                    vec[k] = 0;
                    vec[l] = 1;
                    auto G_l = get_sign_from_ON(vec, l);
                    auto ind_l = get_ind_from_ON(vec, N, ns, fs, buf);
                    cblas_daxpy(Np, T[l*N+k]*G_k*G_l, ci_i, 1, ci_new+ind_l*Np, 1);
                    vec[k] = 1;
                    vec[l] = 0;
                }
            }
        }
        memcpy(ci,ci_new, sizeof(double)*Ns*Np);
//         std::swap(ci_new, ci);
    }
    delete[] buf;
    delete[] vec; 
}

void malmqvist(const int& N, const int& na, const int& nb, double * ci_1, double * ci_2, double * U, double * V){
    // Вспомогательные переменные
    double * U_copy = new double[N*N];
    double * V_copy = new double[N*N];
    
    memcpy(U_copy, U, sizeof(double)*N*N);
    memcpy(V_copy, V, sizeof(double)*N*N);
    
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1)); 
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    int * fa = new int [na * N];
    int * fb = new int [nb * N];
    double * ci_buf = new double[Na*Nb];
    get_vec(na, N, Na, vec_a);
    get_vec(nb, N, Nb, vec_b);
    get_factorials(na, N, fa);
    get_factorials(nb, N, fb);
    // Получаем матрицы для преобразования КВ.
    auto [P_l, T_l] = calc_P_T(N, U_copy);
    auto [P_r, T_r] = calc_P_T(N, V_copy);
    // Меняем КВ коэффициенты по матрице перестановок
    // !!! В конце ci_1 и ci_2 транспонированы
//     printf("S = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_2, 1));
    
    change_ci_with_P(ci_1, ci_2, N, na, Na, vec_a, fa, nb, Nb, vec_b, fb, P_l, P_r);
//     printf("S = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_2, 1));
    // Меняем КВ коэффициенты по матрице T = (1-L)*U^-1 по beta
    change_ci_with_T_one_dim(N, nb, Nb, vec_b, fb, Na, ci_1, ci_buf, T_l);
    change_ci_with_T_one_dim(N, nb, Nb, vec_b, fb, Na, ci_2, ci_buf, T_r);
    // Транспонируем КВ коэффициенты
    transpose(ci_1, ci_2, Nb, Na);
    // Меняем КВ коэффициенты по матрице T = (1-L)*U^-1 по alpha
    change_ci_with_T_one_dim(N, na, Na, vec_a, fa, Nb, ci_1, ci_buf, T_l);
    change_ci_with_T_one_dim(N, na, Na, vec_a, fa, Nb, ci_2, ci_buf, T_r);
//     printf("S_1 = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_1, 1));
//     printf("S_2 = %.10e\n",cblas_ddot(Na*Nb, ci_2, 1, ci_2, 1));
//     printf("S_12 = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_2, 1));
//     getchar();
    // Чистим память
    delete[] ci_buf;
    delete[] vec_a;
    delete[] vec_b;
    delete[] fa; 
    delete[] fb;
    delete[] T_l;
    delete[] T_r;
    delete[] P_l;
    delete[] P_r;
    delete[] V_copy;
    delete[] U_copy;
}

void malmqvist_one_matrix(const int& N, const int& na, const int& nb, double * ci_1, double * U){
    // Вспомогательные переменные
    double * U_copy = new double[N*N];
    memcpy(U_copy, U, sizeof(double)*N*N);
    
    int Na = (int) std::lround(tgammal(N+1) / tgammal(na+1) / tgammal(N-na+1));
    int Nb = (int) std::lround(tgammal(N+1) / tgammal(nb+1) / tgammal(N-nb+1)); 
    int * vec_a = new int [Na*(na+1)];
    int * vec_b = new int [Nb*(nb+1)];
    int * fa = new int [na * N];
    int * fb = new int [nb * N];
    double * ci_buf = new double[Na*Nb];
    get_vec(na, N, Na, vec_a);
    get_vec(nb, N, Nb, vec_b);
    get_factorials(na, N, fa);
    get_factorials(nb, N, fb);
    // Получаем матрицы для преобразования КВ.
    auto [P_l, T_l] = calc_P_T(N, U_copy);
    // Меняем КВ коэффициенты по матрице перестановок
    // !!! В конце ci_1 и ci_2 транспонированы
//     printf("S = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_2, 1));
    
    change_ci_with_P_1(ci_1, N, na, Na, vec_a, fa, nb, Nb, vec_b, fb, P_l);
//     printf("S = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_2, 1));
    // Меняем КВ коэффициенты по матрице T = (1-L)*U^-1 по beta
    change_ci_with_T_one_dim(N, nb, Nb, vec_b, fb, Na, ci_1, ci_buf, T_l);
    // Транспонируем КВ коэффициенты
    transpose_one_matrix(ci_1, Nb, Na);
    // Меняем КВ коэффициенты по матрице T = (1-L)*U^-1 по alpha
    change_ci_with_T_one_dim(N, na, Na, vec_a, fa, Nb, ci_1, ci_buf, T_l);
//     printf("S_1 = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_1, 1));
//     printf("S_2 = %.10e\n",cblas_ddot(Na*Nb, ci_2, 1, ci_2, 1));
//     printf("S_12 = %.10e\n",cblas_ddot(Na*Nb, ci_1, 1, ci_2, 1));
//     getchar();
    // Чистим память
    delete[] ci_buf;
    delete[] vec_a;
    delete[] vec_b;
    delete[] fa; 
    delete[] fb;
    delete[] T_l;
    delete[] P_l;
    delete[] U_copy;
}

