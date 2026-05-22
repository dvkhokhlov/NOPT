#include <stdio.h>


#include "molecule.h"
#include "matr.h"
#include "timer.h"
#include "CI.h"

double ci_t_cutoff=0;
double ci_v_cutoff=0;
double ci_l_cutoff=0;

int init_ci_map(ci_map * m){
#ifdef GOOGLE_D_M
    m[0].set_empty_key(-1);
#endif
    return 1;
}

int init_ci_map_arr(ci_map_arr * m){
#ifdef GOOGLE_D_M
    m[0].set_empty_key(-1);
#endif
    return 1;
}

int init_ci_map_long(ci_map_long * m){
#ifdef GOOGLE_D_M
    m[0].set_empty_key(-1);
#endif
    return 1;
}

int printf_bit(bool b){
    if(b)printf("1");
    else printf("0");
    return 1;
}

int no_eq_occ(ci_key k1,ci_key k2,int n, int n2){
    
//     for(int i=0;i<n;i++)printf("%d",key_to_oc(k1,i));printf("\n");
//     for(int i=0;i<n;i++)printf("%d",key_to_oc(k2,i));printf("\n");
    for(int i=0;i<n;i++)if(k2[i])if(k1[i]){
//         printf("i = %d",i);
        return 0;
    }
    for(int i=0;i<n;i++)if(k2[i+CI_MAX_SPACE])if(k1[i+CI_MAX_SPACE]){
//         printf("i = %d",i);
        return 0;
    }
    
    return 1;
    
}

int printf_ci_key(ci_key k,int n, int p){
    for(int i=0;i<n;i++){
        if(k[i])printf("1");
        else    printf("0");
    }
    for(int i=0;i<p;i++)printf(" ");
    
    printf("|");
    for(int i=CI_MAX_SPACE;i<n+CI_MAX_SPACE;i++){
        if(k[i])printf("1");
        else    printf("0");
    }
    return 0;
}

ci_key ci_link_key(ci_key k1,ci_key k2,int n, int n2){
    
//     printf("ci_link_key (%d,%d):\n",n,n2);
//     getchar();
//     for(int i=0;i<n ;i++)printf("%d",key_to_oc(k1,i));printf("\n");
//     for(int i=0;i<n2;i++)printf("%d",key_to_oc(k2,i));printf("\n");
//     for(int i=0;i<n2;i++)printf("%d",key_to_oc(k1%int(pow(2,n))+k1/int(pow(2,n))*(pow(2,n2))+k2,i));printf("\n");
//     printf("end\n");
//     getchar();
    ci_key kr;//={};
    
    kr.reset();
//     kr.resize(2*n2);
//     printf_ci_key(k1,n,n2-n);printf("\n");
//     printf_ci_key(k2,n2,0);printf("\n");
    for(int i=0   ;i<n   ;i++)kr[i             ] = k1[i]             |k2[i]             ;
    for(int i=n   ;i<n2  ;i++)kr[i             ] =                    k2[i]             ;
    for(int i=0   ;i<n   ;i++)kr[i+CI_MAX_SPACE] = k1[i+CI_MAX_SPACE]|k2[i+CI_MAX_SPACE];
    for(int i=n   ;i<n2  ;i++)kr[i+CI_MAX_SPACE] =                    k2[i+CI_MAX_SPACE];
//     printf("s = %d %d %d\n",k1.size(),k2.size(),kr.size());
//     printf_ci_key(kr,n2,0);printf("\n");
//     getchar();
//     if(res<0){
//         printf("ERROR %d %d\n",n,n2);
//         exit(1);
//     }
    
    
    return kr;//k1%int(pow(2,n))+k1/int(pow(2,n))*(pow(2,n2))+k2;
    
}

int ci_link_sign(ci_key k1,ci_key k2,int n, int n2){
    
//     for(int i=0;i<n ;i++)printf("%d",key_to_oc(k1,i));printf("\n");
//     for(int i=0;i<n2;i++)printf("%d",key_to_oc(k2,i));printf("\n");
//     for(int i=0;i<n2;i++)printf("%d",key_to_oc(ci_link_key(k1,k2,n,n2),i));printf("\n");
    
    int res=0;
    int n_ch=0;
    for(int i=0;i<n;i++){
        if(k2[i])n_ch++;
        if(k1[i])res+=n_ch;
    }
//     printf("resA = %d n_ch = %d\n",res,n_ch);
//     res=0;
    n_ch=0;
//     for(int i=0;i<n ;i++)printf("%d",key_to_oc(k1,i+n ));printf("\n");
//     for(int i=0;i<n2;i++)printf("%d",key_to_oc(k2,i+n2));printf("\n");
//     for(int i=0;i<n2;i++)printf("%d",key_to_oc(ci_link_key(k1,k2,n,n2),i+n2));printf("\n");
    for(int i=0;i<n;i++){
        if(k2[i+CI_MAX_SPACE])n_ch++;
        if(k1[i+CI_MAX_SPACE])res+=n_ch;
    }
    
//     printf("resB = %d n_ch = %d; %d\n",res,n_ch,res%2);
//     getchar();
    if(res%2)return -1;
    return 1;
    
}

ci_key det_key(int * a, int n){
    
    ci_key k/* = {}*/;
    k.reset();
//     for(int i=0; i<z; i++)k.push_back(0);
    for(int i=0; i<n; i++) {
        if(a[i])k.set(i);
//         else    k.push_back(0);
        
    }
    for(int i=0; i<n; i++) {
        if(a[i+n])k.set(i+CI_MAX_SPACE);
//         else    k.push_back(0);
        
    }
    
    return k;

}

ci_key_long det_key_long(int * a, int n){
    
    ci_key_long k/* = {}*/;
    k.reset();
//     for(int i=0; i<z; i++)k.push_back(0);
    for(int i=0; i<n; i++) {
        if(a[i])k.set(i+2*CI_MAX_SPACE);///writing into the 2-a field
//         else    k.push_back(0);
        
    }
    for(int i=0; i<n; i++) {
        if(a[i+n])k.set(i+3*CI_MAX_SPACE);///writing into the 2-b field
//         else    k.push_back(0);
        
    }
    
    return k;

}

// ci_map CI_map(molecule * A, int n_act, int f, int i_S){
//     
//     ci_map m;
//     init_ci_map(&m);
// //     ci_map m1;
//     ci_key key;
//     int n_e=A->n_act_el_alp[f];
//     int n = A->n_csf(0,f,i_S);
// //     double norm=1.0;
// //     double sign=1.0;
//     for(int i_CI=0;i_CI<n;i_CI++){
//         m[det_key(A->csf_occup[0][f][i_S][i_CI].data(),n_act)]=A->csf_coef[0][f][i_S][i_CI];
//     }
//     
//     return m;
//     
// }


double det_F_calc(ci_key key,double *E,int n){
    
    double F=0.0;
//     printf("WARNING:det_F_calc must be checked!!!!!\n");
    for(int i=0;i<n;i++){
        F+=(key[i]+key[i+CI_MAX_SPACE])*E[i];///check it!!!!!!!!!!!!!!!!!
//         printf("F+= (%d+%d)*%e - %e\n",key_to_oc(key,i),key_to_oc(key,i+n),E[i],F);
//         getchar();
    }
    
    return F;
}

int gen_1_el_matr(double * O,ci_map * A,ci_map * B, int n_act, int f){
    
    ci_key id_A;
    ci_key id_T;
    ci_key id_B;
//     printf("WARNING:gen_1_el_matr must be checked!!! possible wrong f value\n");
//     double * c;
    double coef_A;
    double coef_B;
    double sign_A;
    double sign_B;
    ci_map::iterator end = (*B).end();
    ci_map::iterator current;
    for(const auto&d:*A){
        id_A=d.first;

        sign_A=1;
        coef_A=d.second;
        for(int i=0;i<n_act;i++)
            if(id_A[i+f]){
                id_T=id_A;
                id_T.reset(i+f);
                sign_B=1;
                    for(int j=0;j<n_act;j++)
                    if(id_T[j+f]==0){
                        id_B=id_T;
                        id_B.set(j+f);
                        current=(*B).find(id_B);
                        if(current!=end){
                            coef_B=(*current).second;
                            O[i*n_act+j]+=coef_A*coef_B*sign_A*sign_B;
                        }
                    }
                    else
                        sign_B=-sign_B;

                sign_A=-sign_A;
            }
    }
   
    return 0;
}

int gen_1_el_matr_arr(double * O,ci_map_arr * A, int nA, ci_map_arr * B, int nB, int n_act, int f){
    
    ci_key id_A;
    ci_key id_T;
    ci_key id_B;
//     printf("WARNING:gen_1_el_matr must be checked!!! possible wrong f value\n");
//     double * c;
    double * coef_A;
    double * coef_B;
    double sign_A;
    double sign_B;
    ci_map_arr::iterator end = (*B).end();
    ci_map_arr::iterator current;
    
    set_zero_matr(O,nA*nB*n_act*n_act);
    for(const auto&d:*A){
        id_A=d.first;

        sign_A=1;
        coef_A=d.second;
        for(int i=0;i<n_act;i++)
            if(id_A[i+f]){
                id_T=id_A;
                id_T.reset(i+f);
                sign_B=1;
                    for(int j=0;j<n_act;j++)
                    if(id_T[j+f]==0){
                        id_B=id_T;
                        id_B.set(j+f);
                        current=(*B).find(id_B);
                        if(current!=end){
                            coef_B=(*current).second;
                            for(int i_s=0;i_s<nA;i_s++)
                            for(int j_s=0;j_s<nB;j_s++)
                                O[((i_s*nB+j_s)*n_act+i)*n_act+j]+=coef_A[i_s]*coef_B[j_s]*sign_A*sign_B;
                        }
                    }
                    else
                        sign_B=-sign_B;

                sign_A=-sign_A;
            }
    }
   
    return 0;
}

ci_map_long m_transform(ci_map_long C, double * U, int n_i, int n_o, int f, double sign){
    
//     printf("m_transform was not properly changeged USE CI_MAX_SPACE!!!!!!\n");
//     exit(0);
    ///working here!!!!!!!!!!!!!
    
    ci_map_long t;
    init_ci_map_long(&t);
    ci_key_long key;
    ci_key_long tmp_key;
    double c;
    double cU;
    int sign_i;
    int sign_j;
    
    for(const auto&d:C){
        key=d.first;
        c  =d.second;
//         for(int i=0;i<n_o;i++)printf_bit(key[i]);
//         printf(" | ");
//         for(int i=0;i<n_o;i++)printf_bit(key[i+n_o]);
//         printf(" | ");
//         for(int i=0;i<n_i;i++)printf_bit(key[i+2*n_o]);
//         printf(" | ");
//         for(int i=0;i<n_i;i++)printf_bit(key[i+2*n_o+n_i]);
//         printf("\n");
//         getchar();
        sign_i=1;
        for(int i=0;i<n_i;i++)
        if(key[i+(f+2)*CI_MAX_SPACE]){
            sign_j=sign;
            for(int j=0;j<n_o;j++)
            if(!key[j+f*CI_MAX_SPACE]){
                cU=c*U[i*n_o+j];
//                 for(int ii=0;ii<4*n;ii++)printf("%d",key_to_oc(key,ii));
//                 printf("\n");
//                 for(int ii=0;ii<4*n;ii++)printf("%d",key_to_oc(key-pow(2,i+f*n_i+2*n)+pow(2,j+f*n_o),ii));
//                 printf("\n");
//                 getchar();
                if(fabs(cU)>ci_t_cutoff){
                    tmp_key=key;
//                     printf_bit(tmp_key[i+f*n_i+2*n_o]);printf("I\n");
                    tmp_key.reset(i+(f+2)*CI_MAX_SPACE);
//                     printf_bit(tmp_key[i+f*n_i+2*n_o]);printf("R\n");
//                     printf_bit(tmp_key[i+f*n_o]);printf("I\n");
                    tmp_key.set(j+f*CI_MAX_SPACE);
//                     printf_bit(tmp_key[i+f*n_o]);printf("S\n");
                    
//                     printf("tmp_key\n");
//                     for(int i=0;i<n_o;i++)printf_bit(tmp_key[i]);
//                     printf(" | ");
//                     for(int i=0;i<n_o;i++)printf_bit(tmp_key[i+n_o]);
//                     printf(" | ");
//                     for(int i=0;i<n_i;i++)printf_bit(tmp_key[i+2*n_o]);
//                     printf(" | ");
//                     for(int i=0;i<n_i;i++)printf_bit(tmp_key[i+2*n_o+n_i]);
//                     printf("\n");
//                     getchar();
                    t[tmp_key]+=cU*sign_i*sign_j;
                }
            }
            else sign_j=-sign_j;
            sign_i=-sign_i;
        }
    }
    
    return t;
}

ci_key ci_key_l2s(ci_key_long k_i){
    
    ci_key k_o;
    k_o.reset();
    for(int i=0;i<2*CI_MAX_SPACE;i++)
        k_o[i]=k_i[i];
    
    return k_o;
}

ci_map ci_map_l2s(ci_map_long m_i){
    
    printf("l2s in size = %d\n",m_i.size());
    
    double max_el=0;
    
    for(const auto&d:m_i)
        if(fabs(d.second)>max_el) max_el=fabs(d.second);
    
    printf("max = %e\n", max_el);
    
    max_el = max_el * ci_t_cutoff;
    
    ci_map m_o;
    init_ci_map(&m_o);
    for(const auto&d:m_i)
        if(fabs(d.second)>max_el)
        m_o[ci_key_l2s(d.first)]=d.second;
    
//     printf("l2s in size = %d\n",m_o.size());
//     getchar();
    return m_o;
}
    
// ci_map transform_CI(molecule * A, int n_act_i, int n_act_o, int f, int i_S, double * U){
//     
//     ci_map_long m;
//     ci_map_long m1;
//     ci_key_long key;
//     init_ci_map_long(&m);
//     init_ci_map_long(&m1);
//     fprintf(stderr," transformation %d -> %d for state %d\n",n_act_i,n_act_o, i_S);
// //     PrintMatr10(U,n_act_i,n_act_o,1);
// //     getchar();
//     int n_e=A->n_act_el_alp[f];
//     int n = A->n_csf(0,f,i_S);
//     double norm=0.0;
//     double sign=1.0;
//     for(int i_CI=0;i_CI<n;i_CI++){
//         m[det_key_long(A->csf_occup[0][f][i_S][i_CI].data(),n_act_i)]=A->csf_coef[0][f][i_S][i_CI];
//     }
//     for(const auto&d:m){
//         norm+=d.second*d.second;
//     }
// //     printf("norm_input (%d) = %e\n",i_S,norm);
//     norm=1.0;
//     
//     
// 
//     while(true){
//         m1=m_transform(m,U,n_act_i,n_act_o,0,sign);
//         sign=-sign;
//         m.clear();
//         m=m1;
//         m1.clear();
// 
//         norm=norm/n_e;
//         n_e--;
//         fprintf(stderr,"\r %d alpha electrons left; size = %ld        ",n_e,m.size());
//         
//         if(n_e==0)break;
//     }
//     
//     
//     for(const auto&d:m){
//         key=d.first;
//         m[key]=d.second*norm;
//     }
//     
//     n_e=A->n_act_el_bet[f];
//     norm=1.0;
//     sign=1.0;
//     while(true){
//         m1=m_transform(m,U,n_act_i,n_act_o,1,sign);
//  
//         sign=-sign;
//         m.clear();
//         m=m1;
//         m1.clear();
// 
//         norm=norm/n_e;
//         n_e--;
//         fprintf(stderr,"\r %d  beta electrons left; size = %ld        ",n_e,m.size());
// 
//         if(n_e==0)break;
//         
//     }
// 
//     
//     
//     for(const auto&d:m){
//         key=d.first;
//         m[key]=d.second*norm;
//     }
//     
//     norm=0;
//     for(const auto&d:m){
//         norm+=d.second*d.second;
//     }
// //     printf("norm (%d) = %.10e\n",i_S,norm);
// //     getchar();
//     fprintf(stderr,"\r");
// //     getchar();
// //     exit(1);
//     return ci_map_l2s(m);;
//     
// }

int ci_link(ci_map * O, ci_map ci1, int n1, ci_map ci2, int n2){
    
    double cc;
    for(const auto&d1:ci1)
    for(const auto&d2:ci2)
    if(no_eq_occ(d1.first,d2.first,n1,n2)){
        cc=d1.second*d2.second;
        if(fabs(cc)>ci_l_cutoff)
        (*O)[ci_link_key(d1.first,d2.first,n1,n2)]+=cc*ci_link_sign(d1.first,d2.first,n1,n2);
    }
    
    
    return 1;
    
}


int ci_skalar_product(double * O, ci_map_arr *I1, int n1, ci_map_arr *I2, int n2){
    
    set_zero_matr(O,n1*n2);
    
    double * c1;
    double * c2;
    
    ci_map_arr::iterator current;
    ci_map_arr::iterator end = (*I2).end();
    for(const auto&d:*I1){
        c1 = d.second;
        current = (*I2).find(d.first);
        if(current!=end){
            c2 = (*current).second;
            for(int i=0;i<n1;i++)
            for(int j=0;j<n2;j++)
                O[i*n2+j]+=c1[i]*c2[j];
                    
        }
    }
    
    return 1;
    
}

int calc_CI_2el(double * O, ci_map_arr * I1, ci_map_arr * I2, double * V, int n, int ld, int s1, int s2, int f){
    
    ci_key id_I;
//     ci_key id_T;
//     ci_key id_O;
    double * c;
    // fprintf(stderr,"start tensor transformation; size = %ld\n", (*I1).size());
    int count=0;
//     (*O).rehash((*I).size()*3);
    double sign_A1;
    double sign_A2;
    double sign_B1;
    double sign_B2;
    double sign;
    double  * cV;
    double * c2;
    int i_s;
    double c_max;
    cV = new double[s1];
    
    ci_map_arr::iterator current;
    ci_map_arr::iterator end = (*I2).end();
    
    for(const auto&d:(*I1)){
        id_I = d.first;
        c = d.second;
        // fprintf(stderr,"\ri = %d /%ld",count, (*I1).size());
        count++;
        sign_A1=1;
        for(int i=0;i<n-1;i++)if(id_I[i+f]){
//             id_T=id_I;
            id_I.reset(i+f);
            sign_A2=sign_A1;
            for(int j=i+1;j<n;j++)if(id_I[j+f]){
                id_I.reset(j+f);
                sign_B1=1;
                for(int k=0;k<n-1;k++)if(!id_I[k+f]){
//                     id_O=id_T;
                    id_I.set(k+f);
                    sign_B2=sign_B1;
                    for(int l=k+1;l<n;l++)if(!id_I[l+f]){
                        id_I.set(l+f);
                        sign = sign_A1*sign_A2*sign_B1*sign_B2;
                        c_max=0;
//                         printf("c_max = %e\n",c_max);
                        for(i_s=0;i_s<s1;i_s++){

                            cV[i_s]=c[i_s]*(V[((i*ld+k)*ld+j)*ld+l]-V[((i*ld+l)*ld+j)*ld+k])*sign;
                            if(c_max<fabs(cV[i_s]))c_max=fabs(cV[i_s]);
                        }
//                         printf("c_max = %e\n",c_max);
//                         getchar();
                        if(c_max>ci_v_cutoff){
                            current = (*I2).find(id_I);
                            if(current!=end){
                                c2 = (*current).second;
                                for(int i_s=0;i_s<s1;i_s++)
                                for(int j_s=0;j_s<s2;j_s++)
                                    O[i_s*s2+j_s]+=cV[i_s]*c2[j_s];
                                        
                            }
                            
                        }
                        id_I.reset(l+f);
                    }
                    else
                        sign_B2=-sign_B2;
                    id_I.reset(k+f);
                }
                else
                     sign_B1=-sign_B1;
                
                id_I.set(j+f);
                sign_A2=-sign_A2;
            
            }
            id_I.set(i+f);
            sign_A1=-sign_A1;
        }
    }
    // fprintf(stderr,"\n");
//     printf("size = %ld\n", (*O).size());
    
    delete[] cV;
    return 1;
}

int calc_CI_2el_AB(double * O, ci_map_arr * I1, ci_map_arr * I2, double * V, int n, int ld, int s1, int s2){


    ci_key id_I;
//     ci_key id_T;
//     ci_key id_O;
    double * c;
    // fprintf(stderr,"start tensor transformation\n");
    int count=0;
    double sign_A1;
    double sign_A2;
    double sign_B1;
    double sign_B2;
    double sign;
    double  * cV;
    double * c2;
    int i_s;
    double c_max;
    cV = new double[s1];
    ci_map_arr::iterator current;
    ci_map_arr::iterator end = (*I2).end();
    
    
    for(const auto&d:(*I1)){
        id_I = d.first;
        c = d.second;
        // fprintf(stderr,"\ri = %d / %ld",count, (*I1).size());
        count++;
        
        sign_A1=1;
        for(int i=0;i<n;i++)if(id_I[i]){
//             id_T=id_I;
            id_I.reset(i);
            sign_A2=1;
            for(int j=0;j<n;j++)if(id_I[j+CI_MAX_SPACE]){
                id_I.reset(j+CI_MAX_SPACE);
                sign_B1=1;
                for(int k=0;k<n;k++)if(!id_I[k]){
//                     id_O=id_T;
                    id_I.set(k);
                    sign_B2=1;
                    for(int l=0;l<n;l++)if(!id_I[l+CI_MAX_SPACE]){
                        id_I.set(l+CI_MAX_SPACE);
                        sign =sign_A1*sign_A2*sign_B1*sign_B2;
                        c_max=0;
                        for(i_s=0;i_s<s1;i_s++){
                            cV[i_s]=c[i_s]*V[((i*ld+k)*ld+j)*ld+l]*sign;
                            if(c_max<fabs(cV[i_s]))c_max=fabs(cV[i_s]);
                        }
                        
                        if(c_max>ci_v_cutoff){
                            current = (*I2).find(id_I);
                            if(current!=end){
                                c2 = (*current).second;
                                for(int i_s=0;i_s<s1;i_s++)
                                for(int j_s=0;j_s<s2;j_s++)
                                    O[i_s*s2+j_s]+=cV[i_s]*c2[j_s];
                                        
                            }
                        }
                        
                        id_I.reset(l+CI_MAX_SPACE);
                    }
                    else
                        sign_B2=-sign_B2;
                    id_I.reset(k);
                }
                else
                    sign_B1=-sign_B1;
                id_I.set(j+CI_MAX_SPACE);
                sign_A2=-sign_A2;
            }
            id_I.set(i);
            sign_A1=-sign_A1;
        }
    }
    // fprintf(stderr,"\n");
//     printf("size = %ld\n", (*O).size());
    
    delete[] cV;
    return 1;
}

int calc_CI_2el_all(double * O, ci_map_arr * I1, ci_map_arr * I2, double * V, int n, int ld, int s1, int s2){
    
    set_zero_matr(O,s1*s2);
    
    calc_CI_2el(O, I1, I2, V, n, ld, s1, s2, 0);
//     fprintf(stderr,"size after AA %ld\n",(*O).size());
    printf_timer("V tensor transformation AA");
    calc_CI_2el(O, I1, I2, V, n, ld, s1, s2, CI_MAX_SPACE);
//     fprintf(stderr,"size after BB %ld\n",(*O).size());
    printf_timer("V tensor transformation BB");
    calc_CI_2el_AB(O, I1, I2, V, n, ld, s1, s2);
    printf_timer("V tensor transformation AB");
//     fprintf(stderr,"size after AB %ld\n",(*O).size());
    
    return 0;
}


int print_ci_map(ci_map * ci, int n){
    
    for(const auto&wf:*ci){
         for(int i=0;i<n;i++)printf_bit(wf.first[i]);
         printf("|");
         for(int i=0;i<n;i++)printf_bit(wf.first[i+CI_MAX_SPACE]);
         printf("   %e\n", wf.second);
    }
    
    return 1;
}

int print_ci_map_arr(ci_map_arr * ci, int s, int n){
    
    for(const auto&wf:*ci){
         for(int i=0;i<n;i++)printf_bit(wf.first[i]);
         printf("|");
         for(int i=0;i<n;i++)printf_bit(wf.first[i+CI_MAX_SPACE]);
         printf("   %e\n", wf.second[s]);
    }
    
    return 1;
}


int ci_array_link(ci_map * O,ci_map ** I,int f, int * S, int * v){
    
    if(f==0) {
        printf("Wrong number of fragments in ci_array_link\n");
        exit(1);
    }
    if(f==1){
//        fprintf(stderr,"linking is just copy... state %d\n",S[0]);
       O[0]=I[0][S[0]];
    }
    else{
//         fprintf(stderr,"ci_array_link is not written fo n_f>1 ;)\n");
//         exit(1);
        ci_map m;
        ci_map m1;
        init_ci_map(&m);
        init_ci_map(&m1);
        int n;
        m=I[0][S[0]];
        n=v[0];
        for(int j=1;j<f;j++){
            // fprintf(stderr,"linking f[%d]S[%d] and f[%d]S[%d]\n",j-1,S[j-1],j,S[j]);
            ci_link(&m1, m, n,I[j][S[j]],n+v[j]);
            m=m1;
            n+=v[j];
            m1.clear();
        }
        O[0]=m;
        
    }

    return 1;
}

int ci_add(ci_map_arr * O, ci_map * I, int i, int i_max){
    
    double * D;
    for(const auto&wf:*I){
        D=(*O)[wf.first];
        if(D==NULL){
            D = new double [i_max];
            (*O)[wf.first]=D;
            set_zero_matr(D,i_max);
        }
        D[i]=wf.second;
    }
    
    return 0;   
}

int ci_arr_clear(ci_map_arr * O){
    
    double* D;
    
    for(const auto&d:*O){
        D=(*O)[d.first];
        if(D!=NULL)
            delete[] D;
    }
    
    (*O).clear();
    return 0;   
}


int ci_array_link_to_array(ci_map_arr * O,ci_map ** I,int f, int * S, int * v, int i, int i_max){
    
    if(f==0) {
        printf("Wrong number of fragments in ci_array_link\n");
        exit(1);
    }
    if(f==1){
//        fprintf(stderr,"linking is just copy... state %d\n",S[0]);
       ci_add(O,&(I[0][S[0]]),i,i_max);
    }
    else{
//         exit(1);
        ci_map m;
        ci_map m1;
        init_ci_map(&m);
        init_ci_map(&m1);
        int n;
        m=I[0][S[0]];
        n=v[0];
        for(int j=1;j<f;j++){
            // fprintf(stderr,"linking f[%d]S[%d] and f[%d]S[%d]\n",j-1,S[j-1],j,S[j]);
            ci_link(&m1, m, n,I[j][S[j]],n+v[j]);
            m=m1;
            n+=v[j];
            m1.clear();
        }
        ci_add(O,&m,i,i_max);
        
    }

    return 1;
}

double ci_den_sum(ci_map_arr * ci, double E, double * eps, int n){
    
    double R=0;
    for(const auto&det:*ci){
        R+=det.second[0]*det.second[0]/(E-det_F_calc(det.first,eps,n));
    }
//     if(R>0){
//         PrintMatr(eps,1,n,0);
//         for(const auto&det:*ci){
//             printf_ci_key(det.first,n,0);
//             printf("  %e - %e\n",E,det_F_calc(det.first,eps,n));
//         }
//         getchar();
//     }
    return R;
}
double V_AB_22_J_calc(double * V, double * D_A, double * D_B, int a0, int am, int b0, int bm, int ld){
    
    double res=0;
    
   for(int i=a0; i<am; i++)
    for(int j=b0; j<bm; j++)
    for(int k=a0; k<am; k++)
    for(int l=b0; l<bm; l++)
        res += D_A[i*ld+k]*D_B[(j-b0+a0)*ld+l-b0+a0]*V[((i*ld+k)*ld+j)*ld+l];
    
//     printf("res = %e\n",res);
//     getchar();
    
    return res;
    
    
}

double V_AB_22_K_calc(double * V, double * D_A, double * D_B, int a0, int am, int b0, int bm, int ld){
    
    double res=0;
    
    for(int i=a0; i<am; i++)
    for(int j=b0; j<bm; j++)
    for(int k=a0; k<am; k++)
    for(int l=b0; l<bm; l++)
        res += D_A[i*ld+k]*D_B[(j-b0+a0)*ld+l-b0+a0]*V[((i*ld+l)*ld+j)*ld+k];
    
//     printf("res = %e\n",res);
//     getchar();
    
    return res;
    
    
}




ci_map gen_1e_d_a(ci_map_arr * ci,double * e, int n){
    
    ci_map res;
    init_ci_map(&res);
    ci_key id;
//     print_ci_map_arr(ci,0,n);
    
    for(const auto&det:*ci){
        for(int i=0;i<n;i++)if(det.first[i]==0){
            id=det.first;
//             printf_ci_key(id,n,0);printf("\n");
            id.set(i);
//             printf_ci_key(id,n,0);printf("\n");
//             getchar();
            res[id]=det_F_calc(id,e,n);
        }
        for(int i=0;i<n;i++)if(det.first[i+CI_MAX_SPACE]==0){
            id=det.first;
            id.set(i+CI_MAX_SPACE);
            res[id]=det_F_calc(id,e,n);
        }
    }
//     print_ci_map(&res,n);
//     getchar();
    
    return res;
}

ci_map gen_1e_d(ci_map * ci,double * e, int n){
    
    ci_map res;
    init_ci_map(&res);
    ci_key id;
//     print_ci_map_arr(ci,0,n);
    
    for(const auto&det:*ci){
        for(int i=0;i<n;i++)if(det.first[i]==0){
            id=det.first;
//             printf_ci_key(id,n,0);printf("\n");
            id.set(i);
//             printf_ci_key(id,n,0);printf("\n");
//             getchar();
            res[id]=det_F_calc(id,e,n);
        }
        for(int i=0;i<n;i++)if(det.first[i+CI_MAX_SPACE]==0){
            id=det.first;
            id.set(i+CI_MAX_SPACE);
            res[id]=det_F_calc(id,e,n);
        }
    }
//     print_ci_map(&res,n);
//     getchar();
    
    return res;
}

// ci_map gen_2e_d(ci_map_arr * ci,double * e, int n){
//     
//     ci_map res;
//     init_ci_map(&res);
//     ci_key id;
// //     print_ci_map_arr(ci,0,n);
//     
//     for(const auto&det:*ci){
//         for(int i=0  ;i<n;i++)if(det.first[i]==0)
//         for(int j=i+1;j<n;j++)if(det.first[j]==0){
//             id=det.first;
//             id.set(i);
//             id.set(j);
//             res[id]=det_F_calc(id,e,n);
//         }
//         for(int i=0  ;i<n;i++)if(det.first[i]==0)
//         for(int j=0  ;j<n;j++)if(det.first[j+CI_MAX_SPACE]==0){
//             id=det.first;
//             id.set(i);
//             id.set(j+CI_MAX_SPACE);
//             res[id]=det_F_calc(id,e,n);
//         }
//         for(int i=0  ;i<n;i++)if(det.first[i+CI_MAX_SPACE]==0)
//         for(int j=i+1;j<n;j++)if(det.first[j+CI_MAX_SPACE]==0){
//             id=det.first;
//             id.set(i+CI_MAX_SPACE);
//             id.set(j+CI_MAX_SPACE);
//             res[id]=det_F_calc(id,e,n);
//         }
//     }
// //     print_ci_map(&res,n);
// //     getchar();
//     
//     return res;
// }

ci_map gen_1e_e_a(ci_map_arr * ci,double * e, int n){
    
    ci_map res;
    init_ci_map(&res);
    ci_key id;
//     print_ci_map_arr(ci,0,n);
    
    for(const auto&det:*ci){
        for(int i=0;i<n;i++)if(det.first[i]!=0){
            id=det.first;
//             printf_ci_key(id,n,0);printf("\n");
            id.reset(i);
//             printf_ci_key(id,n,0);printf("\n");
//             getchar();
            res[id]=det_F_calc(id,e,n);
        }
        for(int i=0;i<n;i++)if(det.first[i+CI_MAX_SPACE]!=0){
            id=det.first;
            id.reset(i+CI_MAX_SPACE);
            res[id]=det_F_calc(id,e,n);
        }
    }
//     print_ci_map(&res,n);
//     getchar();
    return res;
}

ci_map gen_1e_e(ci_map * ci,double * e, int n){
    
    ci_map res;
    init_ci_map(&res);
    ci_key id;
//     print_ci_map_arr(ci,0,n);
    
    for(const auto&det:*ci){
        for(int i=0;i<n;i++)if(det.first[i]!=0){
            id=det.first;
//             printf_ci_key(id,n,0);printf("\n");
            id.reset(i);
//             printf_ci_key(id,n,0);printf("\n");
//             getchar();
            res[id]=det_F_calc(id,e,n);
        }
        for(int i=0;i<n;i++)if(det.first[i+CI_MAX_SPACE]!=0){
            id=det.first;
            id.reset(i+CI_MAX_SPACE);
            res[id]=det_F_calc(id,e,n);
        }
    }
//     print_ci_map(&res,n);
//     getchar();
    return res;
}

int gen_arr_by_map(ci_map_arr * O, ci_map * ci, int n){

//     ci_map_arr res;
    init_ci_map_arr(O);
    for(const auto&det:*ci){
        (*O)[det.first]=new double[n];
        set_zero_matr((*O)[det.first],n);
    }
    
    return 0;
}

int ci_map_vec_lin_tr(ci_map_arr * ci,double * U, int n){
    
    double * tmp;
    tmp = new double[n];
    
    for(const auto&det:*ci){
        set_zero_matr(tmp,n);
        for(int j=0;j<n;j++)
        for(int i=0;i<n;i++)
            tmp[j]+=det.second[i]*U[j*n+i];
        for(int i=0;i<n;i++)
            det.second[i]=tmp[i];
        
    }
//     for(int i=0;i<n;i++){
//         printf("\n\nCI[%d]\n",i);
//         for(const auto&det:*ci){
//             printf_ci_key(det.first,4,0);
//             printf("\t% .6e\n",det.second[i]);
//         }
//         getchar();
//     }
        
    
    delete[] tmp;
    
    return 0;
}

int ci_vec_lin_tr(double * ci,double * U, int n, int Na, int Nb, int ld){
    
    double * tmp;
    tmp = new double[n];
    
    for(int i_CI=0;i_CI<Na;i_CI++)
    for(int j_CI=0;j_CI<Nb;j_CI++){
        set_zero_matr(tmp,n);
        for(int j=0;j<n;j++)
        for(int i=0;i<n;i++)
            tmp[j]+=ci[(i_CI*Nb+j_CI)*ld+i]*U[j*n+i];
        for(int i=0;i<n;i++)
            ci[(i_CI*Nb+j_CI)*ld+i]=tmp[i];
        
    }
//     for(int i=0;i<n;i++){
//         printf("\n\nCI[%d]\n",i);
//         for(int i_CI=0;i_CI<Na;i_CI++)
//         for(int j_CI=0;j_CI<Nb;j_CI++){
// //             printf_ci_key(det.first,n,0);
//             printf("% .6e\n",ci[(i_CI*Nb+j_CI)*ld+i]);
//         }
//         getchar();
//     }
//         
    
    delete[] tmp;
    
    return 0;
}



int ci_arr_norm(ci_map_arr * ci, int n){
    
    double * norm;
    norm = new double[n];
    set_zero_matr(norm,n);
    
    printf("\n\nnnCI[%d]\n",0);
    for(const auto&det:*ci){
        printf_ci_key(det.first,n,0);
        printf("\t% .6e\n",det.second[0]);
    }
    getchar();
    
    for(const auto&det:*ci){
        for(int i=0;i<n;i++)
            norm[i]+=det.second[i]*det.second[i];
    }
    PrintMatr(norm,n,1,1);
    for(int i=0;i<n;i++)norm[i]=sqrt(norm[i]);
    PrintMatr(norm,n,1,1);
    
    
    for(const auto&det:*ci){
        for(int i=0;i<n;i++)
            det.second[i]=det.second[i]/norm[i];
    }
    printf("\n\nnCI[%d]\n",0);
    for(const auto&det:*ci){
        printf_ci_key(det.first,n,0);
        printf("\t% .6e\n",det.second[0]);
    }
    getchar();
        
    delete[] norm;
    
    return 0;
}
