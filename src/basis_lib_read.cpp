#include "molecule.h"
#include "chem_data.h"
#include "inp_out.h"
#include "pseudo_potential.h"
#include "basis_lib_read.h"
#include "common_vars.h"
#include "defaults.h"
#include "etc.h"
#include "keywords.h"

#define max(a,b)  (((a)<(b))?(b):(a))

int read_l(char A, int line_num){
    if(A=='S')return 0;
    if(A=='P')return 1;
    if(A=='D')return 2;
    if(A=='F')return 3;
    if(A=='G')return 4;
    if(A=='H')return 5;
    if(A=='I')return 6;
    if(A=='K')return 7;
    
//     if(A=='I')return -1;
    return -1;
}

char char_l(int A){
    if(A==0)return 'S';
    if(A==1)return 'P';
    if(A==2)return 'D';
    if(A==3)return 'F';
    if(A==4)return 'G';
    if(A==5)return 'H';
    if(A==6)return 'I';
    if(A==7)return 'K';
    
//     if(A=='I')return -1;
    return 'X';
}


int read_l_from_line(char * A){
    if(strstr(A,"S"))return 0;
    if(strstr(A,"P"))return 1;
    if(strstr(A,"D"))return 2;
    if(strstr(A,"F"))return 3;
    if(strstr(A,"G"))return 4;
    if(strstr(A,"H"))return 5;
    if(strstr(A,"I"))return 6;
    if(strstr(A,"K"))return 7;
    
//     if(A=='I')return -1;
    return -1;
}


std::vector<Shell> basis_lib_read_gbs(molecule * M,  const char * lib_name,
                                      int in_lib, int start, 
                                      std::vector<std::vector<double>> * lib_coef, 
                                      std::vector<int> * shell_center, bool pure,
                                      std::vector<double> * energy,
                                      std::vector<int> * is_core){
    
    std::vector<Shell> basis;
    
    int wrong_atom;
    int line_num;
    int l;
    int l2;//for sp
    int size;
    int skipping=0;
    
    recursive_file  lib_file;
    char * atom_name = new char[3];
    char * line      = new char[BUF_LINE_LENGTH];
    char * full_name = new char[BUF_LINE_LENGTH];
    char * warn_string = new char[BUF_LINE_LENGTH];
    
    
    std::vector<double> alpha;
    std::vector<double> coeff;
    std::vector<double> coeff2;//for sp
    Shell sh;
    
    if(lib_coef    !=nullptr)lib_coef[0].resize(0);
    if(shell_center!=nullptr)shell_center[0].resize(0);
    if(energy      !=nullptr) energy[0].resize(0);
    if(is_core     !=nullptr)is_core[0].resize(0);
    
    if(in_lib==0)sprintf(full_name,"%s\0",lib_name);
    else{
        if(strstr(lib_name,".gbs"))sprintf(full_name,"%s/%s\0",NOPT_LIB,lib_name);
        else                       sprintf(full_name,"%s/%s.gbs\0",NOPT_LIB,lib_name);
    }
    
    lib_file.r_open(full_name);
    if(lib_file.file[0]==NULL){
        printf("ERROR: could not open basis file \n%s\n",full_name);
        exit(1);
    }
    
    lib_file.r_gets(line,BUF_LINE_LENGTH);
    while((strstr(line,"****")==NULL)&&(!lib_file.r_eof())&&(key_word_comp(line,basis_end_kw)==0))lib_file.r_gets(line,BUF_LINE_LENGTH);
    
    if(lib_file.r_eof()||(key_word_comp(line,basis_end_kw)==1)){
        printf("ERROR: wrong basis file format\n");
        printf("       make sure you use psi4 format\n");
        printf("       atoms must begin and end with ****\n");
        exit(0);
    }
    for(int i_a=0;i_a<M->n_atoms;i_a++){
        
        if(is_core!=nullptr)skipping=M->n_ecp_electrons[i_a];//for Huckel_guess only
        
        wrong_atom=1;
        line_num=0;
        lib_file.go_to_begining();
        
        atom_name[0]=M->atom_names[i_a][0];
        atom_name[1]=M->atom_names[i_a][1];
        atom_name[2]=M->atom_names[i_a][2];
        
        sprintf(warn_string,"looking for next atom: keyword - ****, atom  name - \"%s\"",atom_name);
        lib_file.r_smart_fgets(line,&line_num,full_name,warn_string);

        if(M->nucl_charges_full[i_a])
        while(wrong_atom&&(!lib_file.r_eof())&&(key_word_comp(line,basis_end_kw)==0)){
            while(strstr(line,"****")==NULL)lib_file.r_smart_fgets(line,&line_num,full_name,warn_string);
            lib_file.r_smart_fgets(line,&line_num,full_name,"reading atom name");
            if(line[0]==' '){
                printf("ERROR: atom %s was not found in file %s\n",atom_name,full_name);
                printf("       line %d doesn't contain atom name\n",line_num);
                exit(0);
            }
            
            if(line[0]==atom_name[0])
            if(line[1]==atom_name[1]){
                wrong_atom=0;
                lib_file.r_smart_fgets(line,&line_num,full_name, "reading shell type");
                while(strstr(line,"****")==NULL){
//                     printf("%s\n",line);
                    if(energy !=nullptr){
//                         getchar();
                        char * en_line = strstr(line,"E=");
                        if(en_line==nullptr){
                            printf("ERROR: line %d in %s\n", line_num, full_name);
                            printf("       doesn't contain proper shell description\n");
                            printf("       energy is not given\n");
                            printf("       expected something like \"S    3  E=??? CORE\"\n");
                            exit(1);
                        }
                            
                        if(skipping==0) energy[0].push_back(atof(en_line+2));
                    }
                    if(is_core!=nullptr){
                        if     (strstr(line,"VAL" )){if(skipping==0)is_core[0].push_back(0);}
                        else if(strstr(line,"CORE")){if(skipping==0)is_core[0].push_back(1);}
                        else{
                            printf("ERROR: line %d in %s\n", line_num, full_name);
                            printf("       doesn't contain proper shell description\n");
                            printf("       orbital must be CORE or VAL\n");
                            printf("       expected something like \"S    3  E=??? CORE\"\n");
                            exit(1);
                        }
                        
                    }

                    l2=-1;
                    l=read_l(line[0],line_num);
                    
                    if(l==-1){
                        printf("ERROR: line %d in %s\n", line_num, full_name);
                        printf("       doesn't contain shell description\n");
                        printf("       orbital moment is not recognized\n");
                        printf("       expected something like \"S    3\"\n");
                        exit(1);
                    }
                    if(l==0){
                        l2=read_l(line[1],line_num);
                    }
                    size=atoi(line+2);
                    if(size<1){
                        printf("ERROR: line %d in %s\n", line_num, full_name);
                        printf("       contains bad shell description\n");
                        printf("       number of GTO in shell is illegal n=%d\n",size);
                        printf("       expected something like \"S    3\"\n");
                        exit(1);
                    }
                    alpha.resize(size);
                    coeff.resize(size);
                    if(l2==1)
                        coeff2.resize(size);
                    for(int i=0; i<size;i++){
                        lib_file.r_smart_fgets(line,&line_num,full_name, "reading shell parameters");
                        D_change(line);
                        
                        std::stringstream in(line);
                        in >> alpha[i];//=atof(line   );
                        in >> coeff[i];//=atof(line+25);
                        if(l2==1){
                            in >> coeff2[i];//=atof(line+25);
                        }

                    }
                    
                    if(shell_center!=nullptr)if(skipping==0)shell_center[0].push_back(i_a);
                    if(lib_coef    !=nullptr)if(skipping==0)lib_coef[0].push_back(coeff);
                    sh = Shell(alpha,{{l, pure, coeff},},{{M->atom_coord[3*i_a], M->atom_coord[3*i_a+1], M->atom_coord[3*i_a+2]}});
//                         std::cout<<sh<<std::endl;
//                         getchar();
                    if (skipping==0)basis.push_back(sh);
                    else skipping=max(skipping-4*l-2,0);
                    if(l2==1){
                        printf("ERROR: SP-shells were not tested\n");
                        printf("       use separate S ans P shells\n");
                        exit(0);
                        // if(shell_center!=nullptr)if(skipping==0)shell_center[0].push_back(i_a);
                        // if(lib_coef    !=nullptr)if(skipping==0)lib_coef[0].push_back(coeff2);
                        // sh = Shell(alpha,{{l2, pure, coeff2},},{{M->atom_coord[3*i_a], M->atom_coord[3*i_a+1], M->atom_coord[3*i_a+2]}});
                        // if (skipping==0)basis.push_back(sh);
                        // else skipping=max(skipping-4*l2-2,0);
                    
                    }
                    

                    lib_file.r_smart_fgets(line,&line_num,full_name, "reading shell type");
                }
                
            }
            while(strstr(line,"****")==NULL)lib_file.r_smart_fgets(line,&line_num,full_name,warn_string);
        }

    }
    
    
    delete[] line;
    delete[] full_name;
    delete[] atom_name;
    
    lib_file.r_close();
    
    return basis;
    
}

std::vector<Shell> basis_lib_read_exp(molecule * M,  const char * lib_name, 
                                      int in_lib, int start, 
                                      std::vector<std::vector<double>> * lib_coef, 
                                      std::vector<int> * shell_center, bool pure){
    
    std::vector<Shell> basis;
    
    int wrong_atom;
    int line_num;
    int l;
    int l2;//for sp
    int size;
    
    recursive_file  lib_file;
    char * atom_name = new char[3];
    char * line      = new char[BUF_LINE_LENGTH];
    char * full_name = new char[BUF_LINE_LENGTH];
    char * warn_string = new char[BUF_LINE_LENGTH];
    
    
    std::vector<double> alpha;
    std::vector<std::vector<double>> coeff;
    std::vector<double> numbers_from_line;
//     std::vector<double> coeff2;//for sp
    Shell sh;
    
    if(in_lib)sprintf(full_name,"%s/%s\0",NOPT_LIB,lib_name);
    else      sprintf(full_name,"%s\0",lib_name);
    lib_file.r_open(full_name);
    if(lib_file.file[0]==NULL){
        printf("ERROR: could not open basis file \n%s\n",full_name);
        exit(1);
    }
    
    for(int i=0;i<start;i++)lib_file.r_gets(line,BUF_LINE_LENGTH);
    lib_file.r_gets(line,BUF_LINE_LENGTH);
    while((strstr(line,"basis")!=line)&&(!lib_file.r_eof())&&(key_word_comp(line,basis_end_kw)==0))lib_file.r_gets(line,BUF_LINE_LENGTH);
    
    
     if(lib_file.r_eof()||(key_word_comp(line,basis_end_kw)==1)){
        printf("ERROR: wrong basis file format\n");
        printf("       make sure you use .exp format\n");
        printf("       atoms must begin and end with keyword basis in the beginning of the line\n");
        exit(0);
    }
    for(int i_a=0;i_a<M->n_atoms;i_a++){
        
        wrong_atom=1;
        line_num=0;
        lib_file.go_to_begining();
        
        atom_name[0]=M->atom_names[i_a][0];
        atom_name[1]=M->atom_names[i_a][1];
        atom_name[2]=M->atom_names[i_a][2];
        
        sprintf(warn_string,"looking for next atom: keyword - basis, atom  name - \"%s\"",atom_name);
        
        for(int i=0;i<start;i++)lib_file.r_gets(line,BUF_LINE_LENGTH);
        lib_file.r_smart_fgets(line,&line_num,full_name,warn_string);

        if(M->nucl_charges_full[i_a])
        while(wrong_atom&&(!lib_file.r_eof())&&(key_word_comp(line,basis_end_kw)==0)){
            while(strstr(line,"basis")!=line)lib_file.r_smart_fgets(line,&line_num,full_name,warn_string);
            lib_file.r_smart_fgets(line,&line_num,full_name,"reading atom name");
            if(line[0]==' '){
                printf("ERROR: atom %s was not found in file %s\n",atom_name,full_name);
                printf("       line %d doesn't contain atom name\n",line_num);
                exit(0);
            }
            
            if(line[0]==atom_name[0])
            if(line[1]==atom_name[1]){
                wrong_atom=0;
                while(strstr(line,"end")==NULL){
//                     printf("%s\n",line);


                    l2=-1;
                    l=read_l(line[3],line_num);
                    
                    if(l==-1){
                        printf("ERROR: line %d in %s\n", line_num, full_name);
                        printf("       doesn't contain shell description\n");
                        printf("       orbital moment is not recognized\n");
//                         printf("       expected something like \"S    3\"\n");
                        exit(1);
                    }
//                     if(l==0){
//                         l2=read_l(line[1],line_num);
//                     }
                    lib_file.r_smart_fgets(line,&line_num,full_name,"reading shell parameters");
                    lib_file.r_smart_fgets(line,&line_num,full_name,"reading shell parameters");
                    lib_file.r_smart_fgets(line,&line_num,full_name,"reading shell parameters");
                    
                    coeff.resize(0);
                    numbers_from_line.resize(0);
                    
//                     D_change(line);
                    std::stringstream in(line);
                    double read_num;
                    while(in>>read_num)numbers_from_line.push_back(read_num);
                    int n_shells = numbers_from_line.size()-1;
                    
                    coeff.resize(n_shells);
                    
                    alpha.resize(0);
                    alpha.push_back(numbers_from_line[0]);
                    
                    for(int i=0; i<n_shells;i++)coeff[i].resize(0);
                    for(int i=0; i<n_shells;i++)coeff[i].push_back(numbers_from_line[i+1]);
                    
                    while(true){
                        lib_file.r_smart_fgets(line,&line_num,full_name,"reading shell parameters");
//                         D_change(line);
                        std::stringstream in2(line);
                        double read_num;
                        numbers_from_line.resize(0);
                        while(in2>>read_num)numbers_from_line.push_back(read_num);
                        if(numbers_from_line.size()==0)break;
                        if((numbers_from_line.size()-1)!=n_shells){
                            printf("ERROR: while reading shell parameters from %s\n",full_name);
                            printf("       lines %d and %d have different lengths\n",line_num-1, line_num);
                            exit(0);
                        }
                        alpha.push_back(numbers_from_line[0]);
                        for(int i=0; i<n_shells;i++)coeff[i].push_back(numbers_from_line[i+1]);
                    
                    }
                    
                    
                    if(shell_center!=nullptr)for(int i=0; i<n_shells;i++)shell_center[0].push_back(i_a);
                    if(lib_coef    !=nullptr)for(int i=0; i<n_shells;i++)lib_coef[0].push_back(coeff[i]);
                    for(int i=0; i<n_shells;i++){
                        sh = Shell(alpha,{{l, pure, coeff[i]},},{{M->atom_coord[3*i_a], M->atom_coord[3*i_a+1], M->atom_coord[3*i_a+2]}});
//                         std::cout<<sh<<std::endl;
//                         getchar();
                        basis.push_back(sh);
                    }
                    
                    
                }
                
            }
            while(strstr(line,"end")==NULL)lib_file.r_smart_fgets(line,&line_num,full_name,warn_string);
        }

    }
    
    
    delete[] line;
    delete[] full_name;
    delete[] atom_name;
    
    lib_file.r_close();
    
    return basis;
    
}


int ecp_lib_read_gbs (molecule * M, int in_lib, int start, const char * lib_name, const char * suffix){
    
    int type=0;
    if(strstr(suffix,"ECP"))type =1;
    
    
//     int wrong_atom;
    int duplicate;
    int line_num;
//     int l;

//     int size;
    
    recursive_file  lib_file;
    char * atom_name = new char[3];
    char * ecp_name  = new char[7];
    char * line      = new char[BUF_LINE_LENGTH];
    char * full_name = new char[BUF_LINE_LENGTH];
    
    if(in_lib==0)sprintf(full_name,"%s\0",lib_name);
    else{
        if(strstr(lib_name,".gbs"))sprintf(full_name,"%s/%s\0",NOPT_LIB,lib_name);
        else                       sprintf(full_name,"%s/%s.gbs\0",NOPT_LIB,lib_name);
    }
    
    lib_file.r_open(full_name);
    if(lib_file.file[0]==NULL){
        printf("ERROR: could not open basis file \n%s\n",full_name);
        exit(1);
    }
//     fprintf(mol_report,"Reading %s from %s\n", suffix,full_name);
    
    for(int i=0;i<start;i++)lib_file.r_gets(line,BUF_LINE_LENGTH);
    lib_file.r_gets(line,BUF_LINE_LENGTH);
    while((strstr(line,"****")==NULL)&&(!lib_file.r_eof())&&(key_word_comp(line,basis_end_kw)==0))lib_file.r_gets(line,BUF_LINE_LENGTH);
    
    if(lib_file.r_eof()||(key_word_comp(line,basis_end_kw)==1)||(strstr(line,"*****")!=NULL)){
        printf("ERROR: wrong %s file format\n",suffix);
        printf("       make sure you use psi4 format\n");
        printf("       atoms must begin and end with ****\n");
        exit(0);
    }
    
    
    for(int i_a=0;i_a<M->n_atoms;i_a++){
        
        line_num=0;
        lib_file.go_to_begining();
        atom_name[0]=M->atom_names[i_a][0];
        atom_name[1]=M->atom_names[i_a][1];
        atom_name[2]=M->atom_names[i_a][2];
        
//         atom_name_by_charge(M->nucl_charges_full[i_a],atom_name);
        if(atom_name[1]==' ')atom_name[1]='\0';
        if(islower(atom_name[0]))atom_name[0]=toupper(atom_name[0]);
        if(islower(atom_name[1]))atom_name[1]=toupper(atom_name[1]);
        sprintf(ecp_name,"%s-%s\0",atom_name,suffix);
        for(int i=0;i<start;i++)lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
        lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
        duplicate=0;
        while(!lib_file.r_eof()&&(key_word_comp(line,ecp_so_end_kw)==0)){
            if(strstr(line,ecp_name)==line){
                if(duplicate){
                    printf("WARNING: duplicate %s definition \"%s\" at line %d\n",suffix,ecp_name,line_num);
                    break;
                }
                duplicate=1;
                pseudo_potential pp;
                
                int n_ecp_el;
                std::stringstream in(strstr(line,suffix)+3);
                in>>pp.n_pot;
                if(type){
                    in>>pp.n_el;
                    pp.n_pot++;//one more for general term
//                     fprintf(mol_report,"Number of frozen electrons on atom %d is %d\n",i_a,pp.n_el);
                }
                else
                    pp.n_el=0;
                
                pp.atom=i_a;
                pp.set_n_pot();
                
                M->nucl_charges_calc[i_a]-=pp.n_el;
                M->n_ecp_electrons  [i_a]+=pp.n_el;
                M->n_el_calc-=pp.n_el;
                
                for(int i_p=0;i_p<pp.n_pot;i_p++){
                    if(type)pp.l[i_p]=i_p-1;
                    else    pp.l[i_p]=i_p;
                    lib_file.r_smart_fgets(line,&line_num,full_name,"looking for PP");// (keyword - potential)");
                    lib_file.r_smart_fgets(line,&line_num,full_name,"reading PP");
                    int n=atoi(line);
                    if(n<1){
                        printf("ERROR: line %d in %s\n", line_num, full_name);
                        printf("       doesn't contain propper potential description\n");
                        exit(1);
                    }
                    pp.set_pot_n_terms(i_p,n);
                    for(int i_g=0;i_g<n;i_g++){
                        lib_file.r_smart_fgets(line,&line_num,full_name,"reading PP");
                        std::stringstream in(line);
                        in>>pp.n[i_p][i_g];
                        in>>pp.alpha[i_p][i_g];
                        in>>pp.coeff[i_p][i_g];
                    }
                }
//                 pp.print(mol_report);
                if(type)
                    M->PP.push_back(pp);
                else
                    M->SO_PP.push_back(pp);
            }
            lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
        }
    
    }
    
    
//     fprintf(mol_report,"______________________________________\n\n\n");

    delete[] atom_name;
    delete[] ecp_name ;
    delete[] line     ;
    delete[] full_name;
    lib_file.r_close();

    return 0;
}

int ecp_lib_read_exp (molecule * M, int in_lib, int start, const char * lib_name, const char * suffix){
    
    int type=0;
    if(strstr(suffix,"ECP"))type =1;
    
    int duplicate;
    int line_num;
//     int l;

//     int size;
    
    recursive_file  lib_file;
    char * atom_name = new char[3];
    char * ecp_name  = new char[4];
    char * line      = new char[BUF_LINE_LENGTH];
    char * full_name = new char[BUF_LINE_LENGTH];
    
    if(type) sprintf(ecp_name,"%s\0","ecp");
    else     sprintf(ecp_name,"%s\0","so");
    
    if(in_lib)sprintf(full_name,"%s/%s\0",NOPT_LIB,lib_name);
    else      sprintf(full_name,"%s\0",lib_name);
        
    lib_file.r_open(full_name);
    if(lib_file.file[0]==NULL){
        printf("ERROR: could not open basis file \n%s\n",full_name);
        exit(1);
    }
//     fprintf(mol_report,"Reading %s from %s\n", suffix,full_name);
    
    for(int i_a=0;i_a<M->n_atoms;i_a++){
        
        line_num=0;
        lib_file.go_to_begining();
        atom_name[0]=M->atom_names[i_a][0];
        atom_name[1]=M->atom_names[i_a][1];
        atom_name[2]=M->atom_names[i_a][2];
        
//         atom_name_by_charge(M->nucl_charges_full[i_a],atom_name);
        for(int i=0;i<start;i++)lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
        lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
//         duplicate=0;
        while(!lib_file.r_eof()&&(key_word_comp(line,ecp_so_end_kw)==0)){
            if(strstr(line,ecp_name)==line){
                lib_file.r_smart_fgets(line,&line_num,full_name,"reading atom name");
                if(strstr(line,atom_name)==line){
                    
//                     if(duplicate){
//                         printf("WARNING: duplicate %s definition \"%s\" at line %d\n",suffix,ecp_name,line_num);
//                         break;
//                     }
//                     duplicate=1;
                    
                    pseudo_potential pp;
                    
                    int n_ecp_el;
                    if(type){
                        std::stringstream in(strstr(line,"nelec")+5);
                        lib_file.r_smart_fgets(line,&line_num,full_name,"reading number of ECP electrons");
                        in>>pp.n_el;
//                         fprintf(mol_report,"Number of frozen electrons on atom %d is %d\n",i_a,pp.n_el);
                    }
                    else
                        pp.n_el=0;
                    
                    pp.atom=i_a;
                    
                    M->nucl_charges_calc[i_a]-=pp.n_el;
                    M->n_ecp_electrons  [i_a]+=pp.n_el;
                    M->n_el_calc-=pp.n_el;
                    int i_p=0;
                    pp.l.push_back(read_l_from_line(line+3));
                    
                    lib_file.r_smart_fgets(line,&line_num,full_name,"looking for PP");
                    pp.n.resize(0);
                    pp.alpha.resize(0);
                    pp.coeff.resize(0);
                    while(true){
                        
                        std::vector<int>r_n;
                        std::vector<double>r_alpha;
                        std::vector<double>r_coeff;
                        int n_el;
                        double alpha_el;
                        double coeff_el;
                        
                        while(true){
//                             printf("%s\n",line);
//                             getchar();
                            std::stringstream in(line);
                            in>>n_el;
                            in>>alpha_el;
                            in>>coeff_el;
                            r_n.push_back(n_el);
                            r_alpha.push_back(alpha_el);
                            r_coeff.push_back(coeff_el);
                            lib_file.r_smart_fgets(line,&line_num,full_name,"reading PP");
                            if(isalpha(line[0]))break;
                        }
//                         printf("%s\n",line);
//                         getchar();
                        pp.n.push_back(r_n);
                        pp.alpha.push_back(r_alpha);
                        pp.coeff.push_back(r_coeff);
                        i_p++;
                        if(strstr(line,"end"))break;
                        pp.l.push_back(read_l_from_line(line+3));
                        lib_file.r_smart_fgets(line,&line_num,full_name,"reading PP");
                    }
                    pp.n_pot=i_p;
//                     pp.print(mol_report);
//                     getchar();
                    if(type)
                        M->PP.push_back(pp);
                    else
                        M->SO_PP.push_back(pp);
                }
            }
            lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
//             printf("copy: %s\n",line);
        }
    
    }
//     fprintf(mol_report,"______________________________________\n\n\n");

    delete[] atom_name;
    delete[] ecp_name ;
    delete[] line     ;
    delete[] full_name;
    lib_file.r_close();

    return 0;
}

int ecp_lib_read_grpp(molecule * M, int in_lib, int start, const char * lib_name, const char * suffix){
        
    int type=0;
    if(strstr(suffix,"ECP")){
        type =1;
    }
    
    int duplicate;
    int line_num;
//     int l;

//     int size;
    
    recursive_file  lib_file;
    char * atom_name = new char[3];
    char * ecp_name  = new char[4];
    char * line      = new char[BUF_LINE_LENGTH];
    char * full_name = new char[BUF_LINE_LENGTH];
    char * line_shifted;
    
    if(type) sprintf(ecp_name,"%s\0","ecp");
    else     sprintf(ecp_name,"%s\0","so");
    
    
    if(in_lib)sprintf(full_name,"%s/%s\0",NOPT_LIB,lib_name);
    else      sprintf(full_name,"%s\0",lib_name);
    
    lib_file.r_open(full_name);
    if(lib_file.file[0]==NULL){
        printf("ERROR: could not open basis file \n%s\n",full_name);
        exit(1);
    }
//     fprintf(mol_report,"Reading %s from %s\n", suffix,full_name);
    
    for(int i_a=0;i_a<M->n_atoms;i_a++){
        
        line_num=0;
        lib_file.go_to_begining();
        
        atom_name[0]=M->atom_names[i_a][0];
        atom_name[1]=M->atom_names[i_a][1];
        atom_name[2]=M->atom_names[i_a][2];
        
//         atom_name_by_charge(M->nucl_charges_full[i_a],atom_name);
        for(int i=0;i<start;i++)lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
        lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
//         printf("%s\n",line);
//         getchar();
//         duplicate=0;
        while(!lib_file.r_eof()&&(key_word_comp(line,ecp_so_end_kw)==0)){
            if(strstr(line,"************************")!=NULL){
//                 lib_file.r_smart_fgets(line,&line_num,full_name,"reading atom name");
                if(strstr(line,atom_name)==line){
//                     while(isdigit(line[0])==0)line++;
                    int n_noecp_el =atoi(strstr(line,atom_name)+strlen(atom_name));
                    lib_file.r_smart_fgets(line,&line_num,full_name,"reading atom name");
                    
                    for(int i=0;i<3;i++)lib_file.r_smart_fgets(line,&line_num,full_name,"reading atom name");
                    pseudo_potential pp;
//                     std::vector<Shell> out_core;
                    int n_p;
                    int n_sh;
                    int L;
                    int current_shell=0;
                    
                    while(strstr(line,"AREP")==NULL){
                        L = read_l_from_line(line);
                        std::stringstream in(line);
                        in>>n_p;
                        in>>n_sh;
//                         printf("%d %d %d\n",n_p,n_sh,L);
                        std::vector<double> alpha;
                        std::vector<std::vector<double>> coeff;
                        alpha.resize(n_p);
                        coeff.resize(n_sh);
//                         pp.oc_2J.resize(n_sh);
                        line_shifted=line;
                        for(int i=0;i<n_sh;i++){
                            line_shifted=strchr(line_shifted+1,char_l(L));
                            pp.oc_2J.push_back(atoi(line_shifted+1));
//                             printf("%s\n",line_shifted);
//                             printf("%d\n",pp.oc_2J[i]);
//                             getchar();
                        }
                        for(int i=0;i<n_sh;i++)coeff[i].resize(n_p);
                        for(int i=0;i<n_p ;i++){
                            lib_file.r_smart_fgets(line,&line_num,full_name,"reading out_core");
                            std::stringstream in(line);
                            in>>alpha[i];
                            for(int j=0;j<n_sh;j++)in>>coeff[j][i];
                        }
                        
                        for(int i=0;i<n_sh;i++){
                            Shell sh = Shell(alpha,{{L, 1, coeff[i]},},{{M->atom_coord[3*i_a], M->atom_coord[3*i_a+1], M->atom_coord[3*i_a+2]}});
//                             std::cout<<sh<<std::endl;
//                             getchar();
                            pp.oc_lib_shel_coeff.push_back(coeff[i]);
                            pp.out_core.push_back(sh);
                        }
                        
                        lib_file.r_smart_fgets(line,&line_num,full_name,"reading out_core");
//                         printf("%s\n",line);
//                         getchar();
                    }
                    
                    if(type){
                        pp.n_el=M->nucl_charges_calc[i_a]-n_noecp_el;
//                         fprintf(mol_report,"Number of frozen electrons on atom %d is %d\n",i_a,pp.n_el);
                    }
                    else
                        pp.n_el=0;
                    
                    pp.atom=i_a;
                    
                    M->nucl_charges_calc[i_a]-=pp.n_el;
                    M->n_el_calc-=pp.n_el;
                    M->n_ecp_electrons[i_a]+=pp.n_el;
                    int i_p=0;
                    
                    while(strstr(line,"AREP")!=NULL){
//                         printf("AAAAAA:%s\n",line);
                        L = read_l(strstr(line,"AREP")[-2],0);
                        std::stringstream in(line);
                        in>>n_p;
                        in>>n_sh;
//                         printf("%d %d %d\n",n_p,n_sh,L);
                        for(int i=0;i<n_sh;i++)if(pp.out_core[i+current_shell].contr[0].l!=L){
                            printf("ERROR: inconsistent PP and out_core\n");
                            printf("       see line %d and shell %d\n",line_num, i+current_shell);
                            exit(0);
                        }
                        std::vector<double> alpha;
                        std::vector<double> amplitude;
                        std::vector<std::vector<double>> coeff;
                        std::vector<int>n;
                        alpha.resize(n_p);
                        amplitude.resize(n_p);
                        n.resize(n_p);
                        coeff.resize(n_sh);
                        for(int i=0;i<n_sh;i++)coeff[i].resize(n_p);
                        
                        int n_el;
                        double alpha_el;
                        double amp_el;
                        double coeff_el;
                        
//                         printf("%s\n",line);
//                         exit(0);
                        for(int i=0;i<n_p;i++){
                            lib_file.r_smart_fgets(line,&line_num,full_name,"looking for PP");
                            std::stringstream in(line);
                            
                            in>>n_el;
                            n[i]=n_el;
                            
                            in>>alpha_el;
                            alpha[i]=alpha_el;
                            
                            in>>amp_el;
                            if(type==1)amplitude[i]=amp_el;//first number for ecp
                            in>>amp_el;
                            if(type==0)amplitude[i]=amp_el*2.0/(2.0*L+1.0);//second number for so
                            
                            for(int j=0;j<n_sh;j++){
                                in>>coeff_el;
                                coeff[j][i]=coeff_el;
                            }
//                             lib_file.r_smart_fgets(line,&line_num,full_name,"reading PP");
//                             if(isalpha(line[0]))break;
                        }
                        
                        pp.l.push_back(L);
                        pp.n.push_back(n);
                        pp.alpha.push_back(alpha);
                        pp.coeff.push_back(amplitude);
                        
                        for(int j=0;j<n_sh;j++){
                            pp.oc_coeff.push_back(coeff[j]);
                            pp.oc_n.push_back(n);
                            pp.oc_alpha.push_back(alpha);
                        }
                        
                        i_p++;
                        current_shell+=n_sh;
                        lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
                    }
                    if(type)pp.l[pp.l.size()-1]=-1;
                    pp.n_pot=i_p;
                    pp.n_oc=current_shell;
//                     pp.print(mol_report);
                    if(type)
                        M->PP.push_back(pp);
                    else
                        M->SO_PP.push_back(pp);
                }
            }
            if(key_word_comp(line,ecp_so_end_kw)) break;
            lib_file.r_gets(line,BUF_LINE_LENGTH);line_num++;
        }
    
    }
//     fprintf(mol_report,"______________________________________\n\n\n");

    delete[] atom_name;
    delete[] ecp_name ;
    delete[] line     ;
    delete[] full_name;
    lib_file.r_close();

    return 0;
}
