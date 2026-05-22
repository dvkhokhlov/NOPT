
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chem_data.h"
#include "defaults.h"


double covalentR(int n){
       
    // covalent atomic radii, in angstroms
       
    switch(n){
	  case 1 : return 0.38; break;
      case 2 : return 0.32; break;
      case 3 : return 1.34; break;
      case 4 : return 0.9; break;
      case 5 : return 0.82; break;
      case 6 : return 0.77; break;
      case 7 : return 0.75; break;
      case 8 : return 0.73; break;
      case 9 : return 0.71; break;
      case 10 : return 0.69; break;
      case 11 : return 1.54; break;
      case 12 : return 1.3; break;
      case 13 : return 1.18; break;
      case 14 : return 1.11; break;
      case 15 : return 1.06; break;
      case 16 : return 1.02; break;
      case 17 : return 0.99; break;
      case 18 : return 0.97; break;
      case 19 : return 1.96; break;
      case 20 : return 1.74; break;
      case 21 : return 1.44; break;
      case 22 : return 1.36; break;
      case 23 : return 1.25; break;
      case 24 : return 1.27; break;
      case 25 : return 1.39; break;
      case 26 : return 1.25; break;
      case 27 : return 1.26; break;
      case 28 : return 1.21; break;
      case 29 : return 1.38; break;
      case 30 : return 1.31; break;
      case 31 : return 1.26; break;
      case 32 : return 1.22; break;
      case 33 : return 1.19; break;
      case 34 : return 1.16; break;
      case 35 : return 1.14; break;
      case 36 : return 1.1; break;
      case 37 : return 2.11; break;
      case 38 : return 1.92; break;
      case 39 : return 1.62; break;
      case 40 : return 1.48; break;
      case 41 : return 1.37; break;
      case 42 : return 1.45; break;
      case 43 : return 1.56; break;
      case 44 : return 1.26; break;
      case 45 : return 1.35; break;
      case 46 : return 1.31; break;
      case 47 : return 1.53; break;
      case 48 : return 1.48; break;
      case 49 : return 1.44; break;
      case 50 : return 1.41; break;
      case 51 : return 1.38; break;
      case 52 : return 1.35; break;
      case 53 : return 1.33; break;
      case 54 : return 1.3; break;
      case 55 : return 2.25; break;
      case 56 : return 1.98; break;
      case 57 : return 1.69; break;
      case 58 : return 1.65; break;
      case 59 : return 1.65; break;
      case 60 : return 1.65; break;
      case 61 : return 1.65; break;
      case 62 : return 1.65; break;
      case 63 : return 1.65; break;
      case 64 : return 1.65; break;
      case 65 : return 1.65; break;
      case 66 : return 1.65; break;
      case 67 : return 1.65; break;
      case 68 : return 1.65; break;
      case 69 : return 1.65; break;
      case 70 : return 1.65; break;
      case 71 : return 1.6; break;
      case 72 : return 1.5; break;
      case 73 : return 1.38; break;
      case 74 : return 1.46; break;
      case 75 : return 1.59; break;
      case 76 : return 1.28; break;
      case 77 : return 1.37; break;
      case 78 : return 1.28; break;
      case 79 : return 1.44; break;
      case 80 : return 1.49; break;
      case 81 : return 1.48; break;
      case 82 : return 1.47; break;
      case 83 : return 1.46; break;
      case 84 : return 1.46; break;
      case 85 : return 1.46; break;
      case 86 : return 1.45; break;
	  default: printf("warning: no covalent R for element %d\n",n); return 1.5; break;
	};
}

int charge(char * atom, FILE * data){
    
    fseek(data,0,SEEK_SET);
    
    int c=-1;
    char * line = new char[BUF_LINE_LENGTH];
    fgets(line, BUF_LINE_LENGTH, data);
    while(!feof(data)){
        fgets(line, BUF_LINE_LENGTH, data);
        if(line[0]==atom[0])
        if(line[1]==atom[1]){
            c=atoi(line+4);
            break;
        }
    }
    
    delete[] line;
    
    return c;
}

double mass(char * atom, FILE * data){
    
    fseek(data,0,SEEK_SET);
    
    double m=-1.0;
    char * line = new char[BUF_LINE_LENGTH];
    fgets(line, BUF_LINE_LENGTH, data);
    while(!feof(data)){
        fgets(line, BUF_LINE_LENGTH, data);
        if(line[0]==atom[0])
        if(line[1]==atom[1]){
            m=atof(line+12);
            break;
        }
    }
    
    delete[] line;
    
    return m;
}

int ghost(char * atom, FILE * data){
    
    fseek(data,0,SEEK_SET);
    
    int g=0;
    char * line = new char[BUF_LINE_LENGTH];
    fgets(line, BUF_LINE_LENGTH, data);
    while(!feof(data)){
        fgets(line, BUF_LINE_LENGTH, data);
        if(line[0]==atom[0])
        if(line[1]==atom[1]){
            if(strlen(line)>21)
                g=atoi(line+21);
            
            break;
        }
    }
    
    delete[] line;
    
    return g;
}


int atom_name_by_charge(int n,char* name){
       
    // determination of atom name
    name[2]='\0';
    switch(n){
      case 0 :  name[0]='A'; name[1]=' '; return 0; break;
	  case 1 :  name[0]='H'; name[1]=' '; return 0; break;
      case 2 :  name[0]='H'; name[1]='e'; return 0; break;
      case 3 :  name[0]='L'; name[1]='i'; return 0; break;
      case 4 :  name[0]='B'; name[1]='e'; return 0; break;
      case 5 :  name[0]='B'; name[1]=' '; return 0; break;
      case 6 :  name[0]='C'; name[1]=' '; return 0; break;
      case 7 :  name[0]='N'; name[1]=' '; return 0; break;
      case 8  : name[0]='O'; name[1]=' '; return 0; break;
      case 9  : name[0]='F'; name[1]=' '; return 0; break;
      case 10 : name[0]='N'; name[1]=' '; return 0; break;
      case 11 : name[0]='N'; name[1]='a'; return 0; break;
      case 12 : name[0]='M'; name[1]='g'; return 0; break;
      case 13 : name[0]='A'; name[1]='l'; return 0; break;
      case 14 : name[0]='S'; name[1]='i'; return 0; break;
      case 15 : name[0]='P'; name[1]=' '; return 0; break;
      case 16 : name[0]='S'; name[1]=' '; return 0; break;
      case 17 : name[0]='C'; name[1]='l'; return 0; break;
      case 18 : name[0]='A'; name[1]='r'; return 0; break;
      case 19 : name[0]='K'; name[1]=' '; return 0; break;
      case 20 : name[0]='C'; name[1]='a'; return 0; break;
      case 21 : name[0]='S'; name[1]='c'; return 0; break;
      case 22 : name[0]='T'; name[1]='i'; return 0; break;
      case 23 : name[0]='V'; name[1]=' '; return 0; break;
      case 24 : name[0]='C'; name[1]='r'; return 0; break;
      case 25 : name[0]='M'; name[1]='n'; return 0; break;
      case 26 : name[0]='F'; name[1]='e'; return 0; break;
      case 27 : name[0]='C'; name[1]='o'; return 0; break;
      case 28 : name[0]='N'; name[1]='i'; return 0; break;
      case 29 : name[0]='C'; name[1]='u'; return 0; break;
      case 30 : name[0]='Z'; name[1]='n'; return 0; break;
      case 31 : name[0]='G'; name[1]='a'; return 0; break;
      case 32 : name[0]='G'; name[1]='e'; return 0; break;
      case 33 : name[0]='A'; name[1]='s'; return 0; break;
      case 34 : name[0]='S'; name[1]='e'; return 0; break;
      case 35 : name[0]='B'; name[1]='r'; return 0; break;
      case 36 : name[0]='K'; name[1]='r'; return 0; break;
//       case 37 : name[0]='H'; name[1]=' '; return 0; break;
//       case 38 : name[0]='H'; name[1]=' '; return 0; break;
      case 39 : name[0]='Y'; name[1]=' '; return 0; break;
//       case 40 : name[0]='H'; name[1]=' '; return 0; break;
//       case 41 : name[0]='H'; name[1]=' '; return 0; break;
      case 42 : name[0]='M'; name[1]='o'; return 0; break;
//       case 43 : name[0]='H'; name[1]=' '; return 0; break;
//       case 44 : name[0]='H'; name[1]=' '; return 0; break;
//       case 45 : name[0]='H'; name[1]=' '; return 0; break;
//       case 46 : name[0]='H'; name[1]=' '; return 0; break;
//       case 47 : name[0]='H'; name[1]=' '; return 0; break;
//       case 48 : name[0]='H'; name[1]=' '; return 0; break;
      case 49 : name[0]='I'; name[1]='n'; return 0; break;
//       case 50 : name[0]='H'; name[1]=' '; return 0; break;
//       case 51 : name[0]='H'; name[1]=' '; return 0; break;
//       case 52 : name[0]='H'; name[1]=' '; return 0; break;
      case 53 : name[0]='I'; name[1]=' '; return 0; break;
//       case 54 : name[0]='H'; name[1]=' '; return 0; break;
//       case 55 : name[0]='H'; name[1]=' '; return 0; break;
//       case 56 : name[0]='H'; name[1]=' '; return 0; break;
//       case 57 : name[0]='H'; name[1]=' '; return 0; break;
//       case 58 : name[0]='H'; name[1]=' '; return 0; break;
//       case 59 : name[0]='H'; name[1]=' '; return 0; break;
//       case 60 : name[0]='H'; name[1]=' '; return 0; break;
//       case 61 : name[0]='H'; name[1]=' '; return 0; break;
//       case 62 : name[0]='H'; name[1]=' '; return 0; break;
//       case 63 : name[0]='H'; name[1]=' '; return 0; break;
//       case 64 : name[0]='H'; name[1]=' '; return 0; break;
//       case 65 : name[0]='H'; name[1]=' '; return 0; break;
//       case 66 : name[0]='H'; name[1]=' '; return 0; break;
//       case 67 : name[0]='H'; name[1]=' '; return 0; break;
//       case 68 : name[0]='H'; name[1]=' '; return 0; break;
//       case 69 : name[0]='H'; name[1]=' '; return 0; break;
      case 70 : name[0]='Y'; name[1]='b'; return 0; break;
//       case 71 : name[0]='H'; name[1]=' '; return 0; break;
//       case 72 : name[0]='H'; name[1]=' '; return 0; break;
//       case 73 : name[0]='H'; name[1]=' '; return 0; break;
//       case 74 : name[0]='H'; name[1]=' '; return 0; break;
//       case 75 : name[0]='H'; name[1]=' '; return 0; break;
//       case 76 : name[0]='H'; name[1]=' '; return 0; break;
//       case 77 : name[0]='H'; name[1]=' '; return 0; break;
//       case 78 : name[0]='H'; name[1]=' '; return 0; break;
//       case 79 : name[0]='H'; name[1]=' '; return 0; break;
//       case 80 : name[0]='H'; name[1]=' '; return 0; break;
//       case 81 : name[0]='H'; name[1]=' '; return 0; break;
//       case 82 : name[0]='H'; name[1]=' '; return 0; break;
//       case 83 : name[0]='H'; name[1]=' '; return 0; break;
//       case 84 : name[0]='H'; name[1]=' '; return 0; break;
//       case 85 : name[0]='H'; name[1]=' '; return 0; break;
//       case 86 : name[0]='H'; name[1]=' '; return 0; break;
	  default: printf("warning: unknown name for element %d\n",n); return 1; break;
	};
}

double mass_by_charge(int n){
       
    // mass in a.u.
       
    switch(n){
	  case  0 : return  0.000; break;
      case  1 : return  1.008; break;
      case  2 : return  4.002; break;
      case  3 : return  6.941; break;
      case  4 : return  9.012; break;
      case  5 : return 10.081; break;
      case  6 : return 12.011; break;
      case  7 : return 14.007; break;
      case  8 : return 15.999; break;
      case  9 : return 18.998; break;
      case 10 : return 20.179; break;
      case 11 : return 22.990; break;
      case 12 : return 24.305; break;
      case 13 : return 26.982; break;
      case 14 : return 28.086; break;
      case 15 : return 30.974; break;
      case 16 : return 32.066; break;
      case 17 : return 35.453; break;
      case 18 : return 39.948; break;
//       case 19 : return 1.96; break;
//       case 20 : return 1.74; break;
//       case 21 : return 1.44; break;
//       case 22 : return 1.36; break;
//       case 23 : return 1.25; break;
//       case 24 : return 1.27; break;
//       case 25 : return 1.39; break;
//       case 26 : return 1.25; break;
//       case 27 : return 1.26; break;
//       case 28 : return 1.21; break;
//       case 29 : return 1.38; break;
//       case 30 : return 1.31; break;
//       case 31 : return 1.26; break;
//       case 32 : return 1.22; break;
//       case 33 : return 1.19; break;
//       case 34 : return 1.16; break;
//       case 35 : return 1.14; break;
//       case 36 : return 1.1; break;
//       case 37 : return 2.11; break;
//       case 38 : return 1.92; break;
//       case 39 : return 1.62; break;
//       case 40 : return 1.48; break;
//       case 41 : return 1.37; break;
//       case 42 : return 1.45; break;
//       case 43 : return 1.56; break;
//       case 44 : return 1.26; break;
//       case 45 : return 1.35; break;
//       case 46 : return 1.31; break;
//       case 47 : return 1.53; break;
//       case 48 : return 1.48; break;
//       case 49 : return 1.44; break;
//       case 50 : return 1.41; break;
//       case 51 : return 1.38; break;
//       case 52 : return 1.35; break;
//       case 53 : return 1.33; break;
//       case 54 : return 1.3; break;
//       case 55 : return 2.25; break;
//       case 56 : return 1.98; break;
//       case 57 : return 1.69; break;
//       case 58 : return 1.65; break;
//       case 59 : return 1.65; break;
//       case 60 : return 1.65; break;
//       case 61 : return 1.65; break;
//       case 62 : return 1.65; break;
//       case 63 : return 1.65; break;
//       case 64 : return 1.65; break;
//       case 65 : return 1.65; break;
//       case 66 : return 1.65; break;
//       case 67 : return 1.65; break;
//       case 68 : return 1.65; break;
//       case 69 : return 1.65; break;
//       case 70 : return 1.65; break;
//       case 71 : return 1.6; break;
//       case 72 : return 1.5; break;
//       case 73 : return 1.38; break;
//       case 74 : return 1.46; break;
//       case 75 : return 1.59; break;
//       case 76 : return 1.28; break;
//       case 77 : return 1.37; break;
//       case 78 : return 1.28; break;
//       case 79 : return 1.44; break;
//       case 80 : return 1.49; break;
//       case 81 : return 1.48; break;
//       case 82 : return 1.47; break;
//       case 83 : return 1.46; break;
//       case 84 : return 1.46; break;
//       case 85 : return 1.46; break;
//       case 86 : return 1.45; break;
	  default: printf("warning: no atomic mass for element %d\n",n); return 50,0; break;
	};
}
