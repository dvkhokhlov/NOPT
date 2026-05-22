# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "symmetry.h"
# include "version.h"
# include "matr.h"
# include "common_vars.h"

double Unity[9]     = {1, 0, 0,
                       0, 1, 0,
                       0, 0, 1};
                       
double C2x[9]       = {1, 0, 0,
                       0,-1, 0,
                       0, 0,-1};
                    
double C2y[9]       = {-1, 0, 0,
                        0, 1, 0,
                        0, 0,-1};
                    
double C2z[9]       = {-1, 0, 0,
                        0,-1, 0,
                        0, 0, 1};
                        
double sigma_xy[9]  = { 1, 0, 0,
                        0, 1, 0,
                        0, 0,-1};
                        
double sigma_yz[9]  = {-1, 0, 0,
                        0, 1, 0,
                        0, 0, 1};
                        
double sigma_xz[9]  = { 1, 0, 0,
                        0,-1, 0,
                        0, 0, 1};
                        
double Inversion[9] = {-1, 0, 0,
                        0,-1, 0,
                        0, 0,-1};                        

int Cn(double ** Op, int n){
    
    double alpha;
    for(int i=0; i<n-1;i++){
        Op[i]=new double[9];
//         printf("%e\n",M_PI);
        alpha=2*M_PI*(i+1)/n;
        Op[i][0]= cos(alpha);
        Op[i][1]=-sin(alpha);
        Op[i][2]=0.0;
        Op[i][3]=-Op[i][1];
        Op[i][4]= Op[i][0];
        Op[i][5]=0.0;
        Op[i][6]=0.0;
        Op[i][7]=0.0;
        Op[i][8]=1.0;
        
//         PrintMatr(Op[i],3,3,1);
    }
    
    return 0;
}

int Sigma_v(double ** Op, int n){
    
    double alpha;
    for(int i=0; i<n;i++){
        Op[i]=new double[9];
//         printf("%e\n",M_PI);
        alpha=2*M_PI*i/n;
        Op[i][0]= cos(alpha);
        Op[i][1]=-sin(alpha);
        Op[i][2]=0.0;
        Op[i][3]= Op[i][1];
        Op[i][4]=-Op[i][0];
        Op[i][5]=0.0;
        Op[i][6]=0.0;
        Op[i][7]=0.0;
        Op[i][8]=1.0;
        
//         PrintMatr(Op[i],3,3,1);
    }
    
    return 0;
}

int Sigma_d(double ** Op, int n){
    
    double alpha;
    for(int i=0; i<n;i++){
        Op[i]=new double[9];
//         printf("%e\n",M_PI);
        alpha=2*M_PI*(i+0.5)/n;
        Op[i][0]= cos(alpha);
        Op[i][1]= sin(alpha);
        Op[i][2]=0.0;
        Op[i][3]= Op[i][1];
        Op[i][4]=-Op[i][0];
        Op[i][5]=0.0;
        Op[i][6]=0.0;
        Op[i][7]=0.0;
        Op[i][8]=1.0;
        
//         PrintMatr(Op[i],3,3,1);
    }
    
    return 0;
}

int C2(double ** Op, int n){
    
    double alpha;
    for(int i=0; i<n;i++){
        Op[i]=new double[9];
//         printf("%e\n",M_PI);
        alpha=2*M_PI*i/n;
        Op[i][0]= cos(alpha);
        Op[i][1]=-sin(alpha);
        Op[i][2]=0.0;
        Op[i][3]= Op[i][1];
        Op[i][4]=-Op[i][0];
        Op[i][5]=0.0;
        Op[i][6]=0.0;
        Op[i][7]=0.0;
        Op[i][8]=-1.0;
        
//         PrintMatr(Op[i],3,3,1);
    }
    
    return 0;
}


int Sn(double ** Op, int n){
    
    double alpha;
    for(int i=0; i<n-1;i++){
        Op[i]=new double[9];
//         printf("%e\n",M_PI);
        alpha=2*M_PI*(i+1)/n;
        Op[i][0]= cos(alpha);
        Op[i][1]=-sin(alpha);
        Op[i][2]=0.0;
        Op[i][3]=-Op[i][1];
        Op[i][4]= Op[i][0];
        Op[i][5]=0.0;
        Op[i][6]=0.0;
        Op[i][7]=0.0;
        Op[i][8]=-1.0;
        
//         PrintMatr(Op[i],3,3,1);
    }
    
    return 0;
}

                        
symmetry::symmetry(){
    
    group    = nullptr;
    Op_Space = nullptr;
    Op_Orb   = nullptr;
    ch       = nullptr;
    rep_name = nullptr;
    dim      = nullptr;
    P        = nullptr;
    T_D5     = nullptr;
    n_ao_d5  = nullptr;
    
}

int no_group(char *name){
    
    printf("%s version %s does not support %s symmetry\n\n", PROJ_NAME, VERSION, name);
    exit(0);
}

int symmetry::init(char* name){
    
    if(LINEAR){
        name[0]='l';
        name[1]='i';
        name[2]='n';
        T='l';
        n=0;
        t=' ';
        
        n_op = 2;
        n_rep = 11;

        alloc();
        
        sprintf(rep_name[ 0]," S ");
        sprintf(rep_name[ 1]," P-");
        sprintf(rep_name[ 2]," P+");
        sprintf(rep_name[ 3]," D-");
        sprintf(rep_name[ 4]," D+");
        sprintf(rep_name[ 5]," F-");
        sprintf(rep_name[ 6]," F+");
        sprintf(rep_name[ 7]," G-");        
        sprintf(rep_name[ 8]," G+");        
        sprintf(rep_name[ 9]," H-");
        sprintf(rep_name[10]," H+");
        
        group = new char[7];
        
        sprintf(group, "LINEAR\0");
        return 0;
    }
    group = new char[4];
    
    sprintf(group, "%c%c%c\0",name[0],name[1],name[2]);
    
    int l=strlen(name);
    
    if((l!=2)&&(l!=3)){
        printf("wrong length of group name %d; it should be 2 or 3 characters long\n",l);
        no_group(name);
    }
    
    T=name[0];
    
    if ((strcmp(name, "Cs") == 0) || (strcmp(name, "Ci") == 0)) {
        n = 1;
        t = name[1];
    }

    else {
        n = atoi(name + 1);

        if (l == 3)
            t = name[2];
        else
            t = ' ';
    }
    
    // should add printf about conversion of infinity to D2h or C2v 

    // Cs added; mirror plane is xy plane
    if(T == 'C')if(n == 1)if(t == 's') {
        
        n_op = 2;
        n_rep = 2;

        alloc();

        Op_Space[0] = Unity;
        Op_Space[1] = sigma_xy;

        //     for(int i=0;i<n_op;i++)
        //         PrintMatr(Op_Space[i],3,3,0);

        sprintf(rep_name[0],  "A' ");dim[0] =  1.0;ch[0][0] =  1.0;ch[0][1] =  1.0;
        sprintf(rep_name[1], "A\" ");dim[1] =  1.0;ch[1][0] =  1.0;ch[1][1] = -1.0;

        return 0;
    }
    
    // Ci added; TODO give alias S2
    if(T == 'C')if(n == 1)if(t == 'i') {

        n_op = 2;
        n_rep = 2;

        alloc();

        Op_Space[0] = Unity;
        Op_Space[1] = Inversion;

        //     for(int i=0;i<n_op;i++)
        //         PrintMatr(Op_Space[i],3,3,0);

        sprintf(rep_name[0], "Ag ");dim[0] = 1.0;ch[0][0] = 1.0;ch[0][1] =  1.0;
        sprintf(rep_name[1], "Au ");dim[1] = 1.0;ch[1][0] = 1.0;ch[1][1] = -1.0;

        return 0;
    }

    // C2 added; unique direction is z axis
    if(T == 'C')if(n == 2)if(t == ' ') {

        n_op = 2;
        n_rep = 2;

        alloc();

        Op_Space[0] = Unity;
        Op_Space[1] = C2z;

        //     for(int i=0;i<n_op;i++)
        //         PrintMatr(Op_Space[i],3,3,0);

        sprintf(rep_name[0], " A ");dim[0] = 1.0;ch[0][0] = 1.0;ch[0][1] =  1.0;
        sprintf(rep_name[1], " B ");dim[1] = 1.0;ch[1][0] = 1.0;ch[1][1] = -1.0;

        return 0;
    }

    // C2v added; xz and yz mirror planes
    if(T == 'C')if(n == 2)if(t == 'v') {

        n_op = 4;
        n_rep = 4;

        alloc();

        Op_Space[0] = Unity;
        Op_Space[1] = C2z;
        Op_Space[2] = sigma_xz;
        Op_Space[3] = sigma_yz;

        //     for(int i=0;i<n_op;i++)
        //         PrintMatr(Op_Space[i],3,3,0);

        sprintf(rep_name[0], "A1 ");dim[0] = 1.0;ch[0][0] =  1.0;ch[0][1] =  1.0;ch[0][2] =  1.0;ch[0][3] =  1.0;
        sprintf(rep_name[1], "A2 ");dim[1] = 1.0;ch[1][0] =  1.0;ch[1][1] =  1.0;ch[1][2] = -1.0;ch[1][3] = -1.0;
        sprintf(rep_name[2], "B1 ");dim[2] = 1.0;ch[2][0] =  1.0;ch[2][1] = -1.0;ch[2][2] =  1.0;ch[2][3] = -1.0;
        sprintf(rep_name[3], "B2 ");dim[3] = 1.0;ch[3][0] =  1.0;ch[3][1] = -1.0;ch[3][2] = -1.0;ch[3][3] =  1.0;

        return 0;
    }

    // D2 added
    if(T == 'D')if(n == 2)if(t == ' ') {

        n_op = 4;
        n_rep = 4;

        alloc();

        Op_Space[0] = Unity;
        Op_Space[1] = C2z;
        Op_Space[2] = C2y;
        Op_Space[3] = C2x;

        //     for(int i=0;i<n_op;i++)
        //         PrintMatr(Op_Space[i],3,3,0);

        sprintf(rep_name[0], " A ");dim[0] = 1.0;ch[0][0] = 1.0;ch[0][1] =  1.0;ch[0][2] =  1.0;ch[0][3] =  1.0;
        sprintf(rep_name[1], "B1 ");dim[1] = 1.0;ch[1][0] = 1.0;ch[1][1] =  1.0;ch[1][2] = -1.0;ch[1][3] = -1.0;
        sprintf(rep_name[2], "B2 ");dim[2] = 1.0;ch[2][0] = 1.0;ch[2][1] = -1.0;ch[2][2] =  1.0;ch[2][3] = -1.0;
        sprintf(rep_name[3], "B3 ");dim[3] = 1.0;ch[3][0] = 1.0;ch[3][1] = -1.0;ch[3][2] = -1.0;ch[3][3] =  1.0;

        return 0;
    }

    if(T=='D')if(n==2)if(t=='h'){
        
        n_op=8;
        n_rep=8;
        
        alloc();
        
        Op_Space[0]=Unity;
        Op_Space[1]=C2z;
        Op_Space[2]=C2y;
        Op_Space[3]=C2x;
        Op_Space[4]=Inversion;
        Op_Space[5]=sigma_xy;
        Op_Space[6]=sigma_xz;
        Op_Space[7]=sigma_yz;
        
        
    //     for(int i=0;i<n_op;i++)
    //         PrintMatr(Op_Space[i],3,3,0);
        
        sprintf(rep_name[0],"Ag ");dim[0]= 1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;ch[0][4]= 1.0;ch[0][5]= 1.0;ch[0][6]= 1.0;ch[0][7]= 1.0;
        sprintf(rep_name[1],"B1g");dim[1]= 1.0;ch[1][0]= 1.0;ch[1][1]= 1.0;ch[1][2]=-1.0;ch[1][3]=-1.0;ch[1][4]= 1.0;ch[1][5]= 1.0;ch[1][6]=-1.0;ch[1][7]=-1.0;
        sprintf(rep_name[2],"B2g");dim[2]= 1.0;ch[2][0]= 1.0;ch[2][1]=-1.0;ch[2][2]= 1.0;ch[2][3]=-1.0;ch[2][4]= 1.0;ch[2][5]=-1.0;ch[2][6]= 1.0;ch[2][7]=-1.0;
        sprintf(rep_name[3],"B3g");dim[3]= 1.0;ch[3][0]= 1.0;ch[3][1]=-1.0;ch[3][2]=-1.0;ch[3][3]= 1.0;ch[3][4]= 1.0;ch[3][5]=-1.0;ch[3][6]=-1.0;ch[3][7]= 1.0;
        sprintf(rep_name[4],"Au ");dim[4]= 1.0;ch[4][0]= 1.0;ch[4][1]= 1.0;ch[4][2]= 1.0;ch[4][3]= 1.0;ch[4][4]=-1.0;ch[4][5]=-1.0;ch[4][6]=-1.0;ch[4][7]=-1.0;
        sprintf(rep_name[5],"B1u");dim[5]= 1.0;ch[5][0]= 1.0;ch[5][1]= 1.0;ch[5][2]=-1.0;ch[5][3]=-1.0;ch[5][4]=-1.0;ch[5][5]=-1.0;ch[5][6]= 1.0;ch[5][7]= 1.0;
        sprintf(rep_name[6],"B2u");dim[6]= 1.0;ch[6][0]= 1.0;ch[6][1]=-1.0;ch[6][2]= 1.0;ch[6][3]=-1.0;ch[6][4]=-1.0;ch[6][5]= 1.0;ch[6][6]=-1.0;ch[6][7]= 1.0;
        sprintf(rep_name[7],"B3u");dim[7]= 1.0;ch[7][0]= 1.0;ch[7][1]=-1.0;ch[7][2]=-1.0;ch[7][3]= 1.0;ch[7][4]=-1.0;ch[7][5]= 1.0;ch[7][6]= 1.0;ch[7][7]=-1.0;
        
        
        
        return 0;
    }
    
    if(T=='C')if(n==2)if(t=='h'){
        
        n_op=4;
        n_rep=4;
        
        alloc();
        
        Op_Space[0]=Unity;
        Op_Space[1]=C2z;
        Op_Space[2]=Inversion;
        Op_Space[3]=sigma_xy;
    
        
    //     for(int i=0;i<n_op;i++)
    //         PrintMatr(Op_Space[i],3,3,0);
        
        sprintf(rep_name[0],"Ag ");dim[0]= 1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;
        sprintf(rep_name[1],"Bg ");dim[1]= 1.0;ch[1][0]= 1.0;ch[1][1]=-1.0;ch[1][2]= 1.0;ch[1][3]=-1.0;
        sprintf(rep_name[2],"Au ");dim[2]= 1.0;ch[2][0]= 1.0;ch[2][1]= 1.0;ch[2][2]=-1.0;ch[2][3]=-1.0;
        sprintf(rep_name[3],"Bu ");dim[3]= 1.0;ch[3][0]= 1.0;ch[3][1]=-1.0;ch[3][2]=-1.0;ch[3][3]= 1.0;
        
        
        
        return 0;
    }
    
    if(T=='C')if(n==3)if(t==' '){
        
        n_op=3;
        n_rep=2;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,3);
        
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
        sprintf(rep_name[0]," A ");dim[0]=1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;
        sprintf(rep_name[1]," E ");dim[1]=1.0;ch[1][0]= 2.0;ch[1][1]=-1.0;ch[1][2]=-1.0;
        
        
        return 0;
    }
    if(T=='C')if(n==3)if(t=='v'){
        
        n_op=6;
        n_rep=3;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,3);
        Sigma_v(Op_Space+3,3);
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
        sprintf(rep_name[0]," A1");dim[0]=1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;ch[0][4]= 1.0;ch[0][5]= 1.0;
        sprintf(rep_name[1]," A2");dim[1]=1.0;ch[1][0]= 1.0;ch[1][1]= 1.0;ch[1][2]= 1.0;ch[1][3]=-1.0;ch[1][4]=-1.0;ch[1][5]=-1.0;
        sprintf(rep_name[2]," E ");dim[2]=2.0;ch[2][0]= 2.0;ch[2][1]=-1.0;ch[2][2]=-1.0;ch[2][3]= 0.0;ch[2][4]= 0.0;ch[2][5]= 0.0;
        
        
        return 0;
    }
    if(T=='D')if(n==3)if(t==' '){
        
        n_op=6;
        n_rep=3;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,3);
        C2(Op_Space+3,3);
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
        sprintf(rep_name[0]," A1");dim[0]=1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;ch[0][4]= 1.0;ch[0][5]= 1.0;
        sprintf(rep_name[1]," A2");dim[1]=1.0;ch[1][0]= 1.0;ch[1][1]= 1.0;ch[1][2]= 1.0;ch[1][3]=-1.0;ch[1][4]=-1.0;ch[1][5]=-1.0;
        sprintf(rep_name[2]," E ");dim[2]=2.0;ch[2][0]= 2.0;ch[2][1]=-1.0;ch[2][2]=-1.0;ch[2][3]= 0.0;ch[2][4]= 0.0;ch[2][5]= 0.0;
        
        
        return 0;
    }
    if(T=='D')if(n==3)if(t=='h'){
        
        n_op=12;
        n_rep=6;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,3);
        C2(Op_Space+3,3);
        Op_Space[6]=sigma_xy;
        Sn(Op_Space+7,3);
        Sigma_v(Op_Space+9,3);
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
       sprintf(rep_name[0],"A'1") ;dim[0]= 1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;ch[0][4]= 1.0;ch[0][5]= 1.0;ch[0][6]= 1.0;ch[0][7]= 1.0;ch[0][8]= 1.0;ch[0][9]= 1.0;ch[0][10]= 1.0;ch[0][11]= 1.0;
       sprintf(rep_name[1],"A'2") ;dim[1]= 1.0;ch[1][0]= 1.0;ch[1][1]= 1.0;ch[1][2]= 1.0;ch[1][3]=-1.0;ch[1][4]=-1.0;ch[1][5]=-1.0;ch[1][6]= 1.0;ch[1][7]= 1.0;ch[1][8]= 1.0;ch[1][9]=-1.0;ch[1][10]=-1.0;ch[1][11]=-1.0;
       sprintf(rep_name[2],"E' ") ;dim[2]= 2.0;ch[2][0]= 2.0;ch[2][1]=-1.0;ch[2][2]=-1.0;ch[2][3]= 0.0;ch[2][4]= 0.0;ch[2][5]= 0.0;ch[2][6]= 2.0;ch[2][7]=-1.0;ch[2][8]=-1.0;ch[2][9]= 0.0;ch[2][10]= 0.0;ch[2][11]= 0.0;
       sprintf(rep_name[3],"A\"1");dim[3]= 1.0;ch[3][0]= 1.0;ch[3][1]= 1.0;ch[3][2]= 1.0;ch[3][3]= 1.0;ch[3][4]= 1.0;ch[3][5]= 1.0;ch[3][6]=-1.0;ch[3][7]=-1.0;ch[3][8]=-1.0;ch[3][9]=-1.0;ch[3][10]=-1.0;ch[3][11]=-1.0;
       sprintf(rep_name[4],"A\"2");dim[4]= 1.0;ch[4][0]= 1.0;ch[4][1]= 1.0;ch[4][2]= 1.0;ch[4][3]=-1.0;ch[4][4]=-1.0;ch[4][5]=-1.0;ch[4][6]=-1.0;ch[4][7]=-1.0;ch[4][8]=-1.0;ch[4][9]= 1.0;ch[4][10]= 1.0;ch[4][11]= 1.0;
       sprintf(rep_name[5],"E\" ");dim[5]= 2.0;ch[5][0]= 2.0;ch[5][1]=-1.0;ch[5][2]=-1.0;ch[5][3]= 0.0;ch[5][4]= 0.0;ch[5][5]= 0.0;ch[5][6]=-2.0;ch[5][7]= 1.0;ch[5][8]= 1.0;ch[5][9]= 0.0;ch[5][10]= 0.0;ch[5][11]= 0.0;

        
        return 0;
    }

    if(T=='C')if(n==4)if(t=='v'){
        
        n_op=8;
        n_rep=5;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,4);
        Sigma_v(Op_Space+4,2);
        Sigma_d(Op_Space+6,2);
    
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
        sprintf(rep_name[0],"A1 ");dim[0]= 1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;ch[0][4]= 1.0;ch[0][5]= 1.0;ch[0][6]= 1.0;ch[0][7]= 1.0;
        sprintf(rep_name[1],"A2 ");dim[1]= 1.0;ch[1][0]= 1.0;ch[1][1]= 1.0;ch[1][2]= 1.0;ch[1][3]= 1.0;ch[1][4]=-1.0;ch[1][5]=-1.0;ch[1][6]=-1.0;ch[1][7]=-1.0;
        sprintf(rep_name[2],"B1 ");dim[2]= 1.0;ch[2][0]= 1.0;ch[2][1]=-1.0;ch[2][2]= 1.0;ch[2][3]=-1.0;ch[2][4]= 1.0;ch[2][5]= 1.0;ch[2][6]=-1.0;ch[2][7]=-1.0;
        sprintf(rep_name[3],"B2 ");dim[3]= 1.0;ch[3][0]= 1.0;ch[3][1]=-1.0;ch[3][2]= 1.0;ch[3][3]=-1.0;ch[3][4]=-1.0;ch[3][5]=-1.0;ch[3][6]= 1.0;ch[3][7]= 1.0;
        sprintf(rep_name[4],"E  ");dim[4]= 2.0;ch[4][0]= 2.0;ch[4][1]= 0.0;ch[4][2]=-2.0;ch[4][3]= 0.0;ch[4][4]= 0.0;ch[4][5]= 0.0;ch[4][6]= 0.0;ch[4][7]= 0.0;
        
        
        return 0;
    }
    
    
    if(T=='C')if(n==6)if(t=='v'){
        
        n_op=12;
        n_rep=6;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,6);
        Sigma_v(Op_Space+6,3);
        Sigma_d(Op_Space+9,3);
    
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
        sprintf(rep_name[0],"A1 ");dim[0]= 1.0;ch[0][0]= 1.0;ch[0][1]= 1.0;ch[0][2]= 1.0;ch[0][3]= 1.0;ch[0][4]= 1.0;ch[0][5]= 1.0;ch[0][6]= 1.0;ch[0][7]= 1.0;ch[0][8]= 1.0;ch[0][9]= 1.0;ch[0][10]= 1.0;ch[0][11]= 1.0;
        sprintf(rep_name[1],"A2 ");dim[1]= 1.0;ch[1][0]= 1.0;ch[1][1]= 1.0;ch[1][2]= 1.0;ch[1][3]= 1.0;ch[1][4]= 1.0;ch[1][5]= 1.0;ch[1][6]=-1.0;ch[1][7]=-1.0;ch[1][8]=-1.0;ch[1][9]=-1.0;ch[1][10]=-1.0;ch[1][11]=-1.0;
        sprintf(rep_name[2],"B1 ");dim[2]= 1.0;ch[2][0]= 1.0;ch[2][1]=-1.0;ch[2][2]= 1.0;ch[2][3]=-1.0;ch[2][4]= 1.0;ch[2][5]=-1.0;ch[2][6]= 1.0;ch[2][7]= 1.0;ch[2][8]= 1.0;ch[2][9]=-1.0;ch[2][10]=-1.0;ch[2][11]=-1.0;
        sprintf(rep_name[3],"B2 ");dim[3]= 1.0;ch[3][0]= 1.0;ch[3][1]=-1.0;ch[3][2]= 1.0;ch[3][3]=-1.0;ch[3][4]= 1.0;ch[3][5]=-1.0;ch[3][6]=-1.0;ch[3][7]=-1.0;ch[3][8]=-1.0;ch[3][9]= 1.0;ch[3][10]= 1.0;ch[3][11]= 1.0;
        sprintf(rep_name[4],"E1 ");dim[4]= 2.0;ch[4][0]= 2.0;ch[4][1]= 1.0;ch[4][2]=-1.0;ch[4][3]=-2.0;ch[4][4]=-1.0;ch[4][5]= 1.0;ch[4][6]= 0.0;ch[4][7]= 0.0;ch[4][8]= 0.0;ch[4][9]= 0.0;ch[4][10]= 0.0;ch[4][11]= 0.0;
        sprintf(rep_name[5],"E2 ");dim[5]= 2.0;ch[5][0]= 2.0;ch[5][1]=-1.0;ch[5][2]=-1.0;ch[5][3]= 2.0;ch[5][4]=-1.0;ch[5][5]=-1.0;ch[5][6]= 0.0;ch[5][7]= 0.0;ch[5][8]= 0.0;ch[5][9]= 0.0;ch[5][10]= 0.0;ch[5][11]= 0.0;
        
        
        return 0;
    }
    
    if(T=='D')if(n==6)if(t=='h'){
        
        n_op=24;
        n_rep=12;
        
        alloc();
        
        Op_Space[0]=Unity;
        Cn(Op_Space+1,6);
        C2(Op_Space+6,6);
        Op_Space[12]=sigma_xy;
        Sn(Op_Space+13,6);
        Sigma_d(Op_Space+18,3);
        Sigma_v(Op_Space+21,3);
//         for(int i=0; i<n_op;i++){
//             PrintMatr(Op_Space[i],3,3,1);
//         }
        
        sprintf(rep_name[0],"A1g");sprintf(rep_name[1],"A2g");sprintf(rep_name[2],"B1g");sprintf(rep_name[3],"B2g");
        dim[0]    = 1.0;           dim[1]    = 1.0;           dim[2]    = 1.0;           dim[3]    = 1.0;            //dim
        ch [0][0 ]= 1.0;           ch [1][0 ]= 1.0;           ch [2][0 ]= 1.0;           ch [3][0 ]= 1.0;            //E
        ch [0][1 ]= 1.0;           ch [1][1 ]= 1.0;           ch [2][1 ]=-1.0;           ch [3][1 ]=-1.0;            //C6
        ch [0][2 ]= 1.0;           ch [1][2 ]= 1.0;           ch [2][2 ]= 1.0;           ch [3][2 ]= 1.0;            //C3
        ch [0][3 ]= 1.0;           ch [1][3 ]= 1.0;           ch [2][3 ]=-1.0;           ch [3][3 ]=-1.0;            //C2
        ch [0][4 ]= 1.0;           ch [1][4 ]= 1.0;           ch [2][4 ]= 1.0;           ch [3][4 ]= 1.0;            //C3^2
        ch [0][5 ]= 1.0;           ch [1][5 ]= 1.0;           ch [2][5 ]=-1.0;           ch [3][5 ]=-1.0;            //C6^5
        ch [0][6 ]= 1.0;           ch [1][6 ]=-1.0;           ch [2][6 ]= 1.0;           ch [3][6 ]=-1.0;            //C'2 1
        ch [0][7 ]= 1.0;           ch [1][7 ]=-1.0;           ch [2][7 ]=-1.0;           ch [3][7 ]= 1.0;            //C'2 2
        ch [0][8 ]= 1.0;           ch [1][8 ]=-1.0;           ch [2][8 ]= 1.0;           ch [3][8 ]=-1.0;            //C'2 1
        ch [0][9 ]= 1.0;           ch [1][9 ]=-1.0;           ch [2][9 ]=-1.0;           ch [3][9 ]= 1.0;            //C'2 2
        ch [0][10]= 1.0;           ch [1][10]=-1.0;           ch [2][10]= 1.0;           ch [3][10]=-1.0;            //C'2 1
        ch [0][11]= 1.0;           ch [1][11]=-1.0;           ch [2][11]=-1.0;           ch [3][11]= 1.0;            //C'2 2
        ch [0][12]= 1.0;           ch [1][12]= 1.0;           ch [2][12]=-1.0;           ch [3][12]=-1.0;            //sigma_h
        ch [0][13]= 1.0;           ch [1][13]= 1.0;           ch [2][13]= 1.0;           ch [3][13]= 1.0;            //C6  *sigma_h -- S_6
        ch [0][14]= 1.0;           ch [1][14]= 1.0;           ch [2][14]=-1.0;           ch [3][14]=-1.0;            //C3  *sigma_h -- S_3
        ch [0][15]= 1.0;           ch [1][15]= 1.0;           ch [2][15]= 1.0;           ch [3][15]= 1.0;            //C2  *sigma_h -- i
        ch [0][16]= 1.0;           ch [1][16]= 1.0;           ch [2][16]=-1.0;           ch [3][16]=-1.0;            //C3^2*sigma_h -- S_3
        ch [0][17]= 1.0;           ch [1][17]= 1.0;           ch [2][17]= 1.0;           ch [3][17]= 1.0;            //C6^5*sigma_h -- S_6
        ch [0][18]= 1.0;           ch [1][18]=-1.0;           ch [2][18]= 1.0;           ch [3][18]=-1.0;            //sigma_d
        ch [0][19]= 1.0;           ch [1][19]=-1.0;           ch [2][19]= 1.0;           ch [3][19]=-1.0;            //sigma_d
        ch [0][20]= 1.0;           ch [1][20]=-1.0;           ch [2][20]= 1.0;           ch [3][20]=-1.0;            //sigma_d
        ch [0][21]= 1.0;           ch [1][21]=-1.0;           ch [2][21]=-1.0;           ch [3][21]= 1.0;            //sigma_v
        ch [0][22]= 1.0;           ch [1][22]=-1.0;           ch [2][22]=-1.0;           ch [3][22]= 1.0;            //sigma_v
        ch [0][23]= 1.0;           ch [1][23]=-1.0;           ch [2][23]=-1.0;           ch [3][23]= 1.0;            //sigma_v
        
        sprintf(rep_name[4],"E1g");sprintf(rep_name[5],"E2g");
        dim[4]    = 2.0;           dim[5]    = 2.0;           //dim
        ch [4][0 ]= 2.0;           ch [5][0 ]= 2.0;           //E
        ch [4][1 ]= 1.0;           ch [5][1 ]=-1.0;           //C6
        ch [4][2 ]=-1.0;           ch [5][2 ]=-1.0;           //C3
        ch [4][3 ]=-2.0;           ch [5][3 ]= 2.0;           //C2
        ch [4][4 ]=-1.0;           ch [5][4 ]=-1.0;           //C3^2
        ch [4][5 ]= 1.0;           ch [5][5 ]=-1.0;           //C6^5
        ch [4][6 ]= 0.0;           ch [5][6 ]= 0.0;           //C'2 1
        ch [4][7 ]= 0.0;           ch [5][7 ]= 0.0;           //C'2 2
        ch [4][8 ]= 0.0;           ch [5][8 ]= 0.0;           //C'2 1
        ch [4][9 ]= 0.0;           ch [5][9 ]= 0.0;           //C'2 2
        ch [4][10]= 0.0;           ch [5][10]= 0.0;           //C'2 1
        ch [4][11]= 0.0;           ch [5][11]= 0.0;           //C'2 2
        ch [4][12]=-2.0;           ch [5][12]= 2.0;           //sigma_h
        ch [4][13]=-1.0;           ch [5][13]=-1.0;           //C6  *sigma_h -- S_6
        ch [4][14]= 1.0;           ch [5][14]=-1.0;           //C3  *sigma_h -- S_3
        ch [4][15]= 2.0;           ch [5][15]= 2.0;           //C2  *sigma_h -- i
        ch [4][16]= 1.0;           ch [5][16]=-1.0;           //C3^2*sigma_h -- S_3
        ch [4][17]=-1.0;           ch [5][17]=-1.0;           //C6^5*sigma_h -- S_6
        ch [4][18]= 0.0;           ch [5][18]= 0.0;           //sigma_d
        ch [4][19]= 0.0;           ch [5][19]= 0.0;           //sigma_d
        ch [4][20]= 0.0;           ch [5][20]= 0.0;           //sigma_d
        ch [4][21]= 0.0;           ch [5][21]= 0.0;           //sigma_v
        ch [4][22]= 0.0;           ch [5][22]= 0.0;           //sigma_v
        ch [4][23]= 0.0;           ch [5][23]= 0.0;           //sigma_v

        
        
        sprintf(rep_name[6],"A1u");sprintf(rep_name[7],"A2u");sprintf(rep_name[8],"B1u");sprintf(rep_name[9],"B2u");
        dim[6]    = 1.0;           dim[7]    = 1.0;           dim[8]    = 1.0;           dim[9]    = 1.0;            //dim
        ch [6][0 ]= 1.0;           ch [7][0 ]= 1.0;           ch [8][0 ]= 1.0;           ch [9][0 ]= 1.0;            //E
        ch [6][1 ]= 1.0;           ch [7][1 ]= 1.0;           ch [8][1 ]=-1.0;           ch [9][1 ]=-1.0;            //C6
        ch [6][2 ]= 1.0;           ch [7][2 ]= 1.0;           ch [8][2 ]= 1.0;           ch [9][2 ]= 1.0;            //C3
        ch [6][3 ]= 1.0;           ch [7][3 ]= 1.0;           ch [8][3 ]=-1.0;           ch [9][3 ]=-1.0;            //C2
        ch [6][4 ]= 1.0;           ch [7][4 ]= 1.0;           ch [8][4 ]= 1.0;           ch [9][4 ]= 1.0;            //C3^2
        ch [6][5 ]= 1.0;           ch [7][5 ]= 1.0;           ch [8][5 ]=-1.0;           ch [9][5 ]=-1.0;            //C6^5
        ch [6][6 ]= 1.0;           ch [7][6 ]=-1.0;           ch [8][6 ]= 1.0;           ch [9][6 ]=-1.0;            //C'2 1
        ch [6][7 ]= 1.0;           ch [7][7 ]=-1.0;           ch [8][7 ]=-1.0;           ch [9][7 ]= 1.0;            //C'2 2
        ch [6][8 ]= 1.0;           ch [7][8 ]=-1.0;           ch [8][8 ]= 1.0;           ch [9][8 ]=-1.0;            //C'2 1
        ch [6][9 ]= 1.0;           ch [7][9 ]=-1.0;           ch [8][9 ]=-1.0;           ch [9][9 ]= 1.0;            //C'2 2
        ch [6][10]= 1.0;           ch [7][10]=-1.0;           ch [8][10]= 1.0;           ch [9][10]=-1.0;            //C'2 1
        ch [6][11]= 1.0;           ch [7][11]=-1.0;           ch [8][11]=-1.0;           ch [9][11]= 1.0;            //C'2 2
        ch [6][12]=-1.0;           ch [7][12]=-1.0;           ch [8][12]= 1.0;           ch [9][12]= 1.0;            //sigma_h
        ch [6][13]=-1.0;           ch [7][13]=-1.0;           ch [8][13]=-1.0;           ch [9][13]=-1.0;            //C6  *sigma_h -- S_6
        ch [6][14]=-1.0;           ch [7][14]=-1.0;           ch [8][14]= 1.0;           ch [9][14]= 1.0;            //C3  *sigma_h -- S_3
        ch [6][15]=-1.0;           ch [7][15]=-1.0;           ch [8][15]=-1.0;           ch [9][15]=-1.0;            //C2  *sigma_h -- i
        ch [6][16]=-1.0;           ch [7][16]=-1.0;           ch [8][16]= 1.0;           ch [9][16]= 1.0;            //C3^2*sigma_h -- S_3
        ch [6][17]=-1.0;           ch [7][17]=-1.0;           ch [8][17]=-1.0;           ch [9][17]=-1.0;            //C6^5*sigma_h -- S_6
        ch [6][18]=-1.0;           ch [7][18]= 1.0;           ch [8][18]=-1.0;           ch [9][18]= 1.0;            //sigma_d
        ch [6][19]=-1.0;           ch [7][19]= 1.0;           ch [8][19]=-1.0;           ch [9][19]= 1.0;            //sigma_d
        ch [6][20]=-1.0;           ch [7][20]= 1.0;           ch [8][20]=-1.0;           ch [9][20]= 1.0;            //sigma_d
        ch [6][21]=-1.0;           ch [7][21]= 1.0;           ch [8][21]= 1.0;           ch [9][21]=-1.0;            //sigma_v
        ch [6][22]=-1.0;           ch [7][22]= 1.0;           ch [8][22]= 1.0;           ch [9][22]=-1.0;            //sigma_v
        ch [6][23]=-1.0;           ch [7][23]= 1.0;           ch [8][23]= 1.0;           ch [9][23]=-1.0;            //sigma_v
        
        sprintf(rep_name[10],"E1u");sprintf(rep_name[11],"E2u");
        dim[10]    = 2.0;           dim[11]    = 2.0;           //dim
        ch [10][0 ]= 2.0;           ch [11][0 ]= 2.0;           //E
        ch [10][1 ]= 1.0;           ch [11][1 ]=-1.0;           //C6
        ch [10][2 ]=-1.0;           ch [11][2 ]=-1.0;           //C3
        ch [10][3 ]=-2.0;           ch [11][3 ]= 2.0;           //C2
        ch [10][4 ]=-1.0;           ch [11][4 ]=-1.0;           //C3^2
        ch [10][5 ]= 1.0;           ch [11][5 ]=-1.0;           //C6^5
        ch [10][6 ]= 0.0;           ch [11][6 ]= 0.0;           //C'2 1
        ch [10][7 ]= 0.0;           ch [11][7 ]= 0.0;           //C'2 2
        ch [10][8 ]= 0.0;           ch [11][8 ]= 0.0;           //C'2 1
        ch [10][9 ]= 0.0;           ch [11][9 ]= 0.0;           //C'2 2
        ch [10][10]= 0.0;           ch [11][10]= 0.0;           //C'2 1
        ch [10][11]= 0.0;           ch [11][11]= 0.0;           //C'2 2
        ch [10][12]= 2.0;           ch [11][12]=-2.0;           //sigma_h
        ch [10][13]= 1.0;           ch [11][13]= 1.0;           //C6  *sigma_h -- S_6
        ch [10][14]=-1.0;           ch [11][14]= 1.0;           //C3  *sigma_h -- S_3
        ch [10][15]=-2.0;           ch [11][15]=-2.0;           //C2  *sigma_h -- i
        ch [10][16]=-1.0;           ch [11][16]= 1.0;           //C3^2*sigma_h -- S_3
        ch [10][17]= 1.0;           ch [11][17]= 1.0;           //C6^5*sigma_h -- S_6
        ch [10][18]=-0.0;           ch [11][18]=-0.0;           //sigma_d
        ch [10][19]=-0.0;           ch [11][19]=-0.0;           //sigma_d
        ch [10][20]=-0.0;           ch [11][20]=-0.0;           //sigma_d
        ch [10][21]=-0.0;           ch [11][21]=-0.0;           //sigma_v
        ch [10][22]=-0.0;           ch [11][22]=-0.0;           //sigma_v
        ch [10][23]=-0.0;           ch [11][23]=-0.0;           //sigma_v
 
        return 0;
    }

    
    printf("error: group %c%d%c is not found\n",T,n,t);
    
    exit(0);
    
    return 0;
    
}


int symmetry::alloc(){
    
    if(LINEAR==0)Op_Space = new double*[n_op ];
    if(LINEAR==0)Op_Orb   = new double*[n_op ];
    if(LINEAR==1)T_D5     = new double*[n_rep];
    if(LINEAR==1)n_ao_d5  = new int    [n_rep];
    
    ch       = new double*[n_rep];
    dim      = new double [n_rep];
    P        = new double*[n_rep];for(int i=0;i<n_rep;i++)P[i]=nullptr;
    rep_name = new char  *[n_rep];
    for(int i=0;i<n_rep;i++)
        rep_name[i]=new char[4];
    
    for(int i=0;i<n_rep;i++){
        ch[i]=new double[n_op];
    }
    
    return 0;
}

int symmetry::print_rep(FILE * out, int n){
    
    if(n==-1)
        fprintf(out,"       A   ");
    else
        fprintf(out,"      %s  ",rep_name[n]);
    
    return 0;
    
}


symmetry::~symmetry(){
    
    if(group    != nullptr) delete[]group;
    
    if(Op_Space != nullptr) delete[]Op_Space;
    
    if(Op_Orb   != nullptr){for(int i=0;i<n_rep;i++)delete[]Op_Orb  [i];delete[]Op_Orb  ;}
    if(ch       != nullptr){for(int i=0;i<n_rep;i++)delete[]ch      [i];delete[]ch      ;}
    if(rep_name != nullptr){for(int i=0;i<n_rep;i++)delete[]rep_name[i];delete[]rep_name;}
    if(T_D5     != nullptr){for(int i=0;i<n_rep;i++)delete[]T_D5    [i];delete[]T_D5    ;}
    if(dim      != nullptr) delete[]dim;
    if(n_ao_d5  != nullptr) delete[]n_ao_d5;
    
    if(P        != nullptr){
        for(int i=0;i<n_rep;i++)if(P[i]!= nullptr)delete[]P[i];
        delete[]P;
    }
}


int OxR(double * OR, double * O, double *R){
    
    OR[0]=O[0]*R[0]+O[1]*R[1]+O[2]*R[2];
    OR[1]=O[3]*R[0]+O[4]*R[1]+O[5]*R[2];
    OR[2]=O[6]*R[0]+O[7]*R[1]+O[8]*R[2];
    
    return 0;
}

int at_search(double * R, double * L, int n, double eps){
    
    for(int i=0;i<n;i++){
        if(fabs(R[0]-L[3*i+0])<eps)
        if(fabs(R[1]-L[3*i+1])<eps)
        if(fabs(R[2]-L[3*i+2])<eps)
            return i;
    }
    
    return -1;
    
}
