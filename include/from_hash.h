//#include <google/dense_hash_map>
#include <bitset>
#include "CI.h"
#include "molecule.h"



inline void get_factorials(const int& n, const int& N, int * f){
    for (int i=0; i<n; i++)
        for (int j=1; j<N; j++)
            f [ i*N + j ] = (int) std::lround(tgammal(j+1) / tgammal(i+1) / tgammal(j-i+1));
}

// Генерирует индекс для bitset ON вектора. Нужно передать вспомогательные расчитанные факториалы f и 
// int * буферный массив размера (n+1).
inline int get_ind_from_ON(int * vec, const int N, const int n, int * f, int * buf){
    int count = 1;
    for (int j=0; j<N; j++)
        if (vec[j]){
            buf[count]=j+1;
            count+=1;
        }
    int p = 0;
    for (int k=1; k<=n; k++)
        for (int j=buf[k-1]+1; j<buf[k]; j++){
            p+=f[ (n-k)*N + N-j ];
        }
    return p;
}

inline int get_index_from_array(int * occ, const int& n, const int& N, int * f, int * buf, const int& increment){
    int count = 1;
    for (int j=0; j<N; j++)
        if (occ[j+increment]){
            buf[count]=j+1;
            count+=1;
        }
    int p = 0;
    for (int k=1; k<=n; k++)
        for (int j=buf[k-1]+1; j<buf[k]; j++){
            p+=f[ (n-k)*N + N-j ];
        }
    return p;
}

inline int get_sign_from_ON(int * vec, const int& k){
    int p = 0;
    for (int j=0; j<k; j++) p += vec[j];
    return (p%2) ? -1 : 1;
}

// Поиск детерминантов с занятой k орбиталью и все однократные возбуждения по k-ой орбитали
inline void get_I_I_friends(const int& n, const int& N, const int& size, int * f, const int& k, int * v, int * Ia, int * Ia_friends){
    int * buf = new int[n+1];
    int * vec = new int[N];
    buf[0] = 0;
    int count = 0;
    for (int p=0; p<size; p++){
        // Обнуляем вектор
        for (int i=0; i<N; i++) vec[i]=0;
        // Заполняем ON вектор
        for (int i=1; i<n+1; i++) { vec[v[p*(n+1)+i]-1] = 1;}
        // Проходимся по занятым
        if (vec[k]) {
            auto G_k = get_sign_from_ON(vec, k);
            Ia[count] = p;
            int internal_count = 0;
            for (int j=0; j<N; j++){
                if (!vec[j]){
                    vec[k] = 0;
                    vec[j] = 1;
                    Ia_friends[count*(N-n)*3 + internal_count*3 + 0] = get_ind_from_ON(vec, N, n, f, buf);
                    Ia_friends[count*(N-n)*3 + internal_count*3 + 1] = j;
                    Ia_friends[count*(N-n)*3 + internal_count*3 + 2] = get_sign_from_ON(vec, j) * G_k;
                    internal_count += 1;
                    vec[k] = 1;
                    vec[j] = 0;
                }
            }
            count+=1;
        }
    }
    delete[] buf;
    delete[] vec;
}



inline int get_index(const ci_key& key, const int& n, const int& N, int * f, int * buf, const int& increment){
    int count = 1;
    for (int j=0; j<N; j++)
        if (key[j+increment]){
            buf[count]=j+1;
            count+=1;
        }
    int p = 0;
    for (int k=1; k<=n; k++)
        for (int j=buf[k-1]+1; j<buf[k]; j++){
            p+=f[ (n-k)*N + N-j ];
        }
    return p;
}

// void get_I_I_friends(const int& n, const int& N, const int& size, int * f, const int& k, int * v, int * Ia, int * Ia_friends);

void LU(const int& size, double * A, int * P);

void InvU(const int& size, double * U);

void change_ci_with_P_one_dim(double * ci, const int& N, const int& ns, const int& Ns, int * vec_s, int * fs, const int& Np, int * P);

void change_ci_with_T_one_dim(const int& N, const int& ns, const int& Ns, int * vec_s, int * fs, const int& Np, double *& ci, double *& ci_new, double * T);

void get_vec(const int& n, const int& N, const int& size, int * v);

void ci_from_hash(const int& N, const int& na, const int& nb, ci_map_arr& ci_hash, const int& i_state, double *& ci);

void hash_from_ci(const int& N, int n_empty, const int& na, const int& nb, ci_map & ci_hash, const int& i_state, double *& ci);

void occ_from_ci(const int& N, const int& na, const int& nb, std::vector<std::vector<int>> * occ);

void occ_from_ci_ba(const int& N, const int& na, const int& nb, std::vector<std::vector<int>> * occ);

void ci_from_mol(const int& N, const int& na, const int& nb, molecule * M, int s, int f, const int& i_state, double *& ci);

double calc_2e_ints(double * H, const int& N, const int& na, const int& nb, double * ci_1, double * ci_2, double * G,
                    const bool& singlet, const int n_s1, const int ld1, const int n_s2, const int ld2);

int gen_1_el_DMA_aldet(double * O,double * ci_A, double * ci_B, int nA, int nB, int n_act, int na, int nb);

int gen_1_el_DMB_aldet(double * O,double * ci_A, double * ci_B, int nA, int nB, int n_act, int na, int nb);

//Функция для преобразования КВ коэффициентов при пересчете актиных молекулярных орбиталей
////      Принимает параметры пространства N, na, nb, КВ коэффициенты в виде матриц,
////      и две матрицы преобразования левых и правых активных молекулярных орбиталей.
////      !!!! Не VT, а уже V. Для транспонирования есть функция ниже

void malmqvist(const int& N, const int& na, const int& nb, double * ci_1, double * ci_2, double * Left, double * Right);

void malmqvist_one_matrix(const int& N, const int& na, const int& nb, double * ci_1, double * U);

void transpose_one_matrix(double * A, const int& rows, const int& columns);
