//C headers
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <vector>
# include "common_vars.h"
# include "grabbers.h"
# include "defaults.h"
# include "inp_out.h"
# include "keywords.h"


recursive_file::recursive_file(){
    
    file.resize(0);
    level = 0;
}

int recursive_file::r_open(char * name){
    
    level =0;
    file.resize(1);
    file[0] = fopen(name,"r");
    
    return 0;
}

int recursive_file::go_to_begining(){
    
    level=0;
    file.resize(1);
    return fseek(file[0],0,SEEK_SET);
    
}


int recursive_file::r_gets(char * out, int l){
    
    fgets(out, l, file[level]);
    
    while(feof(file[level])){
        if(level==0)return 0;
        level--;
        fgets(out, l, file[level]);
    }
        
    char * new_file;
    FILE * F;
    
    while(key_word_comp(out, include_kw)){
//         printf("%s\n",out);
        kw_to_s(&new_file, out, include_kw);
        F=fopen(new_file,"r");
        if(F!=nullptr){
            level++;
            file.resize(level+1);
            file[level]=F;
            fgets(out, l, file[level]);
            while(feof(file[level])){
                level--;
                fgets(out, l, file[level]);
            }
        }
        else{
            printf("ERROR: can not find included file %s.\n", new_file);
            exit(0);
        }
            
//         printf("%s\n",out);
//         getchar();
    }
    
    
    return 1;
    
}

int recursive_file::r_smart_fgets(char * l, int * n, char* fn, const char * job){
    
    r_gets(l,BUF_LINE_LENGTH);
    n[0]++;
    if(r_eof()){
       printf("ERROR: unexpected end of file %s\n",fn);
       printf("       at line %d\n",n[0]);
       printf("       while %s\n", job);
       exit(0);
    }
    
    return 0;
}

int recursive_file::r_eof(){
    if(level)return 0;
    if(feof(file[0]))return 1;
    
    return 0;
}

int recursive_file::r_close(){
    
    for(int i=0;i<file.size();i++)
        fclose(file[i]);
    
    file.resize(0);
    level = 0;
    
    return 0;
}

recursive_file::~recursive_file(){
    
    r_close();
    
}

int smart_fgets(char * l, int * n, FILE* f, char* fn, const char * job){
    
    fgets(l,BUF_LINE_LENGTH,f);
    n[0]++;
    if(feof(f)){
       printf("ERROR: unexpected end of file %s\n",fn);
       printf("       at line %d\n",n[0]);
       printf("       while %s\n", job);
       exit(0);
    }
    
    return 0;
}


int print_file_content(char * inp){
    
    recursive_file f;
    
    f.r_open(inp);
    
    if(f.file[0]==NULL){
        fprintf(out_stream,"ERROR: could not open %s\n",inp);
        exit(1);
    }
    
    
    fprintf(out_stream,"\nInput file:\n\n");
    char s[BUF_LINE_LENGTH];
    
    int i=0;
    f.r_gets(s,BUF_LINE_LENGTH);
    while(!f.r_eof()){
        printf("%4d: %s",i,s);
        f.r_gets(s,BUF_LINE_LENGTH);i++;
    }
    fprintf(out_stream,"_______________________________________________________________________\n");
    
    
    return 0;
}

int read_block_size(const char * B){
    
    FILE * in;
    char fname[50];
    sprintf(fname,"%s%s%s",B,B,".out\0");
    in = fopen(fname,"r");
    if(in==NULL){
        fprintf(out_stream,"couldn't open data for %s and %s\n",B,B);
        return -1;
    }
    
    char buf[255];
    fgets(buf,255,in);
    int result2 = 0;
    
    while(!feof(in)){
        
        if(strstr(buf,"S[")){
            fgets(buf,255,in);
            fgets(buf,255,in);
            fgets(buf,255,in);
            if(strstr(buf,"H[")){
                fgets(buf,255,in);
                fgets(buf,255,in);
                if(strstr(buf,"D[")){
                    result2++;
                    fgets(buf,255,in);
                }
            }
        }
        
        fgets(buf,255,in);
    
    }
    
    int result = sqrt(result2);
    
    if(result2-result*result){
        fprintf(out_stream,"ERROR: number of elements in %s is not square (n_el = %d)\n",fname, result2);
        return -1;
    }
    
    return result;
}

int read_HS_data(const char * B, const char * K, double * H, double * S, double * x, double * y, double * z, int ny, int nx, int ld){
    
    FILE * in;
    char fname[50];
    sprintf(fname,"%s%s%s",B,K,".out\0");
    in = fopen(fname,"r");
    int di,dj,ii,jj;
    int ii_old=-1;
    int jj_old=-1;
    int line_was_switched=0;
    di = ld-nx+1;
    dj = 1;
    if(in==NULL){
        int tmp = ny;
        ny=nx;
        nx=tmp;
        dj=ld;
        di=1-ld*(nx-1);
        sprintf(fname,"%s%s%s",K,B,".out\0");
        in = fopen(fname,"r");
        if(in==NULL){
            fprintf(out_stream,"couldn't open data for %s and %s\n",B,K);
            return 1;
        }
    }
    
//     fprintf(out_stream," %d %d %d %d\n", nx,ny,di,dj);
    
//     getchar();
    char buf[255];
    fgets(buf,255,in);
    int i=0;
    int j=0;
    while(!feof(in)){
        
        if(strstr(buf,"S[")){
            sscanf(&(buf[2]),"%ld,%ld] = %lf",&ii,&jj,S);
//             fprintf(out_stream,"start reading %d %d (%d %d) ((%d %d)) %e from %s\n",ii,jj,ii,jj,ii_old,jj_old, S[0],buf);
            fgets(buf,255,in);
            fgets(buf,255,in);
            fgets(buf,255,in);
            if(strstr(buf,"H[")){
                sscanf(&(buf[2]),"%ld,%ld] = %lf",&ii,&jj,H);
                fgets(buf,255,in);
                fgets(buf,255,in);
                if(strstr(buf,"D[")){
                    sscanf(&(buf[2]),"%ld,%ld] = {%lf, %lf,%lf",&ii,&jj,x,y,z);
                    //control of order of indexes
                    if(line_was_switched){
                        if(jj>jj_old)if(ii<ii_old){
                            fprintf(out_stream,"Wrong(type1) element index in file %s. May be name is reversed.\n",fname,nx,ny);
                            return 1;
                        }
                    }
                    else{
                        if(jj<jj_old){
                            fprintf(out_stream,"Wrong(type2) element index in file %s. May be name is reversed.\n",fname,nx,ny);
                            return 1;
                        }
                    }
                    //prepare next step
                    j++;
                    if(j==nx){
                        line_was_switched=1;
                        j=0;
                        i++;
                        S+=di;
                        H+=di;
                        x+=di;
                        y+=di;
                        z+=di;
                        
                        
                    }
                    else{
                        line_was_switched=0;
                        S+=dj;
                        H+=dj;
                        x+=dj;
                        y+=dj;
                        z+=dj;
                    }
                    
                    if(i==ny)if(j==1){
                            fprintf(out_stream,"block overflow while reading file %s for block %dx%d\n",fname,nx,ny);
                            return 1;
                        }
                    fgets(buf,255,in);
                    
                }
                else{
                    fprintf(out_stream,"D[%d,%d] not found in file %s\n",ii,jj,fname);
                    return 1;
                }
            }
            else{
                fprintf(out_stream,"H[%d,%d] not found in file %s\n",ii,jj,fname);
                return 1;
            }
            ii_old=ii;
            jj_old=jj;
        }
        
        fgets(buf,255,in);
    
    }
    
    if((i==ny)&&(j==0)){
        return 0;
    }
    
    fprintf(out_stream,"insuffisient data while reading file %s for block %dx%d. Element [%d,%d] not found.\n",fname,nx,ny,i,j);
    return 1;
    
//     fprintf(out_stream,"Data for <%s|%s>:\n",B,K);
//     fprintf(out_stream,"H = %.10f\n",H[0]);
//     fprintf(out_stream,"S = %.10f\n",S[0]);
//     fprintf(out_stream,"D = {%.10f, %.10f, %.10f}\n",x[0],y[0],z[0]);

}

int read_H_data(char * B, double * H, int nx){
    
    FILE * in;
    char fname[50];
    sprintf(fname,"%s%s%s",B,"v_xmc",".out\0");
    in = fopen(fname,"r");
    int di,dj,ii,jj;
    int ii_old=-1;
    int jj_old=-1;
    int line_was_switched=0;
    di = 1;
    dj = 1;

    
//     fprintf(out_stream," %d %d %d %d\n", nx,ny,di,dj);
    
//     getchar();
    char buf[255];
    fgets(buf,255,in);
    int i=0;
    int j=0;
    while(!feof(in)){
        
        if(strstr(buf,"H[")){
            sscanf(&(buf[2]),"%ld,%ld] = %lf",&ii,&jj,H);
//             fprintf(out_stream,"start reading %d %d (%d %d) ((%d %d)) %e from %s\n",ii,jj,ii,jj,ii_old,jj_old, S[0],buf);
            //control of order of indexes
            if(line_was_switched){
                if(jj>jj_old)if(ii<ii_old){
                    fprintf(out_stream,"Wrong(type1) element index in file %s. Check consistency of doCI and PT data.\n",fname);
                    return 1;
                }
            }
            else{
                if(jj<jj_old){
                    fprintf(out_stream,"Wrong(type2) element index in file %s. Check consistency of doCI and PT data.\n",fname);
                    return 1;
                }
            }
            //prepare next step
            j++;
            if(j==nx){
                line_was_switched=1;
                j=0;
                i++;
            }
            
            H++;
            
            ii_old=ii;
            jj_old=jj;
        }
        
        fgets(buf,255,in);
    
    }
    
    if((i==nx)&&(j==0)){
        return 0;
    }
    
    fprintf(out_stream,"insuffisient data while reading file %s for block %dx%d. Element [%d,%d] not found.\n",fname,nx,nx,i,j);
    return 1;
    
//     fprintf(out_stream,"Data for <%s|%s>:\n",B,K);
//     fprintf(out_stream,"H = %.10f\n",H[0]);
//     fprintf(out_stream,"S = %.10f\n",S[0]);
//     fprintf(out_stream,"D = {%.10f, %.10f, %.10f}\n",x[0],y[0],z[0]);

} 

int read_gv_eigenvalue(double *C, char * ev_name, int * dim, int w, int s){
    
    FILE * f;
    char str[1023];
    char state_n[12];
    
    f = fopen(ev_name,"r");
    if(f==NULL){
        fprintf(out_stream,"Unable to open source file %s\n",ev_name);
        exit(1);
    }
    
    fgets(str,1023,f);
    
    while((strstr(str,"Eigenvectors:")==0)&&(!feof(f))) fgets(str,1023,f);
    
    if (feof(f)){
        fprintf(out_stream,"Didn't find \"Eigenvectors:\" in %s\n",ev_name);
    }
    
    sprintf(state_n,"State %d:",s);

    while((strstr(str,state_n)==0)&&(!feof(f))) fgets(str,1023,f);
    
    if (feof(f)){
        fprintf(out_stream,"Didn't find \"%s\" in %s\n",state_n,ev_name);
    }
    
    //TODO: make protection!!!
    int i_s=0;
    for(int i=0;i<w;i++){
        fgets(str,1023,f);
//         fprintf(out_stream,"%s",str);
//         fprintf(out_stream,"%s",strstr(str,"|")+1);
        for(int j=0;j<dim[i];j++){
//             fprintf(out_stream,"%s",strstr(str,"|")+14*j+1);
            sscanf(strstr(str,"|")+14*j+1,"%lf",C+i_s);
            fprintf(out_stream,"% .4f ",C[i_s]);
            i_s++;
        }
        fprintf(out_stream,"\n");
    }
    
    return 1;
}

int ReadMatr10(double * M, int n, int m, char * f_name, const char * keyword){
    
    FILE * inp;
    inp = fopen(f_name,"r");
    if(inp==NULL){
        fprintf(out_stream,"WARNING: input file %s was not found by ReadMatr10\n",f_name);
        return 1;
    }
    char line[255];
    fgets(line,255,inp);
    while(strstr(line,keyword)==NULL){
        fgets(line,255,inp);
        if(feof(inp)){
            fprintf(out_stream,"%s keyword id not found in %s\n", keyword, f_name);
            exit(1);
        }
    }
    fgets(line,255,inp);
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
           M[i*m+j]=atof(line+19*j);
        }
        fgets(line,255,inp);
    }
    
    fclose(inp);
    
    return 0;
    
}

int ReadDipole(double * x, double * y, double * z, int n, char * f_name, const char * keyword){
    
    FILE * inp;
    inp = fopen(f_name,"r");
    if(inp==NULL){
        fprintf(out_stream,"ERROR:input file was not found\n");
        fprintf(out_stream,"      by function ReadDipole\n");
        exit(1);
    }
    
    char line[255];
    fgets(line,255,inp);
    while(strstr(line,keyword)==NULL){
        fgets(line,255,inp);
        if(feof(inp)){
            fprintf(out_stream,"%s keyword id not foubnd in %s\n", keyword, f_name);
            exit(1);
        }
    }
    
    fgets(line,255,inp);
    fgets(line,255,inp);
    
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            x[i*n+j]=atof(line+14);
            y[i*n+j]=atof(line+28);
            z[i*n+j]=atof(line+42);
            fgets(line,255,inp);
//             fprintf(out_stream,"%e %e %e\n",x[i*n+j],y[i*n+j],z[i*n+j]);
//             getchar();
        }
        fgets(line,255,inp);
    }
    
    
    fclose(inp);
    
    return 0;
    
}

int PrintDipoleTable(FILE * stream, int dim){
    
    for(int i=0;i<dim;i++)
        fprintf(stream,"                                          |             |");
    fprintf(stream,"\n");
    
    for(int i=0;i<dim;i++)
        fprintf(stream,"                d[%3d]                    |     |d|     |",i,i);
    
    fprintf(stream,"\n");
    
    for(int i=0;i<dim;i++)
        fprintf(stream,"__________________________________________|_____________|");
    fprintf(stream,"\n");
    fflush(stream);
    
    return 0;
}

int PrintEnergyTable(FILE * stream, int dim){
    
    fprintf(stream,"                  | ");
    for(int i=1;i<dim;i++)
        fprintf(stream,"                  |           | ");
    fprintf(stream,"\n");
    
    fprintf(stream,"  E[  0]          | ");
    for(int i=1;i<dim;i++)
        fprintf(stream,"E[%3d]            | dE[%3d]   | ",i,i);
    
    fprintf(stream,"\n");
    fprintf(stream,"__________________| ");
    for(int i=1;i<dim;i++)
        fprintf(stream,"__________________|___________| ");
    fprintf(stream,"\n");
    fflush(stream);
    
    return 0;
}


int PrintDipole(double * x, double * y, double * z, int n){
    
    for(int i=0;i<n;i++){
        fprintf(out_stream,"\n");
        for(int j=0;j<n;j++){
            fprintf(out_stream,"d[%3d,%3d] = {% 5e,% 5e,% 5e} |d| = %5e\n",i,j,
                                                               x[i*n+j],
                                                               y[i*n+j],
                                                               z[i*n+j],
                                                          sqrt(x[i*n+j]*x[i*n+j]+        
                                                               y[i*n+j]*y[i*n+j]+
                                                               z[i*n+j]*z[i*n+j]));
        }
    }
    
    if(grab_D){
        
        if(n>grabber_n_states){
            nopt_printf("ERROR: can not copy dipole of size %dx%d\n",n,n);
            nopt_printf("       by the grabber of size %dx%d\n",grabber_n_states,grabber_n_states);
            exit(0);
        }
        memcpy(Dx_grabbed,x,n*n*sizeof(double));
        memcpy(Dy_grabbed,y,n*n*sizeof(double));
        memcpy(Dz_grabbed,z,n*n*sizeof(double));
    }
    
    if(d_scan!=nullptr){
        for(int i=0;i<n;i++)
            fprintf(d_scan,"% 5e,% 5e,% 5e | %.5e |",x[i*n+i],
                                                     y[i*n+i],
                                                     z[i*n+i],
                                                sqrt(x[i*n+i]*x[i*n+i]+        
                                                     y[i*n+i]*y[i*n+i]+
                                                     z[i*n+i]*z[i*n+i]));
        fprintf(d_scan,"\n");
        fflush(d_scan);
    }
    
    fprintf(out_stream,"\n\n");
    
    return 1;
    
}

int PrintDispersion(double * xx, double * yy, double * zz, int n){
    
    for(int i=0;i<n;i++){
        fprintf(out_stream,"\n");
        fprintf(out_stream,"State %d:\n",i);
        fprintf(out_stream,"<x2> = % 5e\n",xx[i*(n+1)]);
        fprintf(out_stream,"<y2> = % 5e\n",yy[i*(n+1)]);
        fprintf(out_stream,"<z2> = % 5e\n",zz[i*(n+1)]);
        
    }
    
    fprintf(out_stream,"\n\n");
    
    return 1;
    
}

int PrintQuadrupole(double * xx, double * yy, double * zz,
                    double * xy, double * xz, double * yz, int n){
    
    for(int i=0;i<n;i++){
        fprintf(out_stream,"\n");
        for(int j=0;j<n;j++){
            fprintf(out_stream,"Q[%3d,%3d] = {% 5e,% 5e,% 5e \n",i,j,
                                                     xx[i*n+j],
                                                     xy[i*n+j],
                                                     xz[i*n+j]);
            fprintf(out_stream,"              % 5e,% 5e,% 5e \n",xy[i*n+j],
                                                     yy[i*n+j],
                                                     yz[i*n+j]);
            fprintf(out_stream,"              % 5e,% 5e,% 5e}\n",xz[i*n+j],
                                                     yz[i*n+j],
                                                     zz[i*n+j]);
        }
    }
    
    fprintf(out_stream,"\n\n");
    
    return 1;
    
}

int PrintEnergy(double * ev, int dim, int scan){
    
    fprintf(out_stream,"E[  0] = %.10f\n",ev[0]);
    for(int i=1;i<dim;i++)
        fprintf(out_stream,"E[%3d] = %.10f | %.4f eV\n",i,ev[i],27.2114*(ev[i]-ev[0]));
//         fprintf(out_stream,"E[%3d] = %.10f | w[%3d] = %7.2f nm | %.4f eV\n",i,ev[i],i,45.5634/(ev[i]-ev[0]),27.2114*(ev[i]-ev[0]));
    
    if(scan)if(e_scan!=nullptr){
        fprintf(e_scan,"%.10f H | ",ev[0]);
        for(int i=1;i<dim;i++)
            fprintf(e_scan,"%.10f H | %.4f eV | ",ev[i],27.2114*(ev[i]-ev[0]));
        fprintf(e_scan,"\n");
        fflush(e_scan);
    }
    
    if(grab_E==1){
        E_grabbed.resize(dim);
        for(int i=0;i<dim;i++) E_grabbed[i]=ev[i];
    }
  
    

    return 0;
}

int PrintEnergy_derivative_corrected(double * E, double * D, int dim, int change){
    
    fprintf(out_stream,"____________________________________\n");
    fprintf(out_stream," State| dE/dx        | rel.err. (%%) |\n");
    fprintf(out_stream,"______|______________|______________|\n");
    for(int i=0;i<dim;i++)
        fprintf(out_stream,"  %3d | % .5e |   % 6.2f     |\n",i,D[i*(dim+1)],(D[i*(dim+1)]-D[0])*100);
    fprintf(out_stream,"______|______________|______________|\n\n");
    
//     double * E_pr = new double[dim];
//     
//     for(int i=0;i<dim;i++)
//         E_pr[i]=E[0]+(E[i]-E[0])*(1+D[i*(dim+1)]-D[0]);
//     
//     PrintEnergy(E_pr,dim);
//     
//     if(change)
//         for(int i=0;i<dim;i++)
//             E[i]=E_pr[i];
//     
//     delete[] E_pr;
//     
    
    return 0;
}

int PrintDer(double * D, int n_s, int n_d){
    
    fprintf(out_stream,"______");
    for(int i_d=0;i_d<n_d;i_d++)
    fprintf(out_stream,"_______________");
    fprintf(out_stream,"\n");
    fprintf(out_stream," State|");
    for(int i_d=0;i_d<n_d;i_d++)
    fprintf(out_stream," dE/dx(%3d)   |",i_d);
    fprintf(out_stream,"\n");
    fprintf(out_stream,"______|");
    for(int i_d=0;i_d<n_d;i_d++)
    fprintf(out_stream,"______________|");
    fprintf(out_stream,"\n");
    for(int i=0;i<n_s;i++){
        double c =0;
        fprintf(out_stream,"  %3d |",i);
        for(int i_d=0;i_d<n_d;i_d++)
        fprintf(out_stream," % .5e |",D[i*n_d+i_d]);
        for(int i_d=0;i_d<n_d;i_d++)
            c+=D[i*n_d+i_d];
        fprintf(out_stream," % .5e |",c);
        fprintf(out_stream,"\n");
    }
    fprintf(out_stream,"______|");
    for(int i_d=0;i_d<n_d;i_d++)
    fprintf(out_stream,"______________|");
    fprintf(out_stream,"\n");
    
    return 0;
}


int PrintRowVec(double * V, int n1, int n2){
    for(int i=0;i<n1;i++){
        fprintf(out_stream,"v[%3d] =",i);
        for(int j=0;j<n2;j++)fprintf(out_stream," % .10e ",V[i*n2+j]);
        fprintf(out_stream,"\n");
    }
    return 0;
}


std::vector<double> read_E_XMC(const char* B, const char * E){
    
    std::vector<double> result;
    result.resize(0);
    FILE * in;
    char fname[50];
    sprintf(fname,"%s%s%s",B,E,"\0");
    in = fopen(fname,"r");
    if(in==NULL){
        fprintf(out_stream,"couldn't open %s\n",fname);
        exit(1);;
    }
    
    
    double a;
    char * buf=new char[255];
    char * e;
    fgets(buf,255,in);
    while(!feof(in)){
        e=strstr(buf,"E(MP2)=");
        fprintf(out_stream,"%s\n",buf);
//         fprintf(out_stream,"%c\n",e[0]);
        if(e!=NULL){
            a=atof(e+7);
            result.push_back(a);
        }
        
        fgets(buf,255,in);
    }
//     fprintf(out_stream,"%d\n",result.size());
//     for(int i=0;i<result.size();i++)
//         fprintf(out_stream,"E%s[%d]=%e\n",E,i,result[i]);
    
    
    return result;
    
}

int print_CI_results(FILE * target, int i, int j, double S, double V_nuc, double H_1el, double H_2el, double d_x, double d_y, double d_z){
    
    fprintf(target,"S[%d,%d] = %.10f\n",i,j,S);
    fprintf(target,"H1[%d,%d] = %.10f\n",i,j,H_1el+V_nuc*S);
    fprintf(target,"H1[%d,%d]/S[%d,%d] = %.10f\n",i,j,i,j,H_1el/S+V_nuc);
    fprintf(target,"H[%d,%d] = %.10f\n",i,j,H_1el+H_2el+V_nuc*S);
    fprintf(target,"H2[%d,%d]/S[%d,%d] = %.10f\n",i,j,i,j,(H_1el+H_2el)/S+V_nuc);
    fprintf(target,"D[%d,%d] = {%.10f, %.10f, %.10f}\n",i,j,d_x,d_y,d_z);
    fprintf(target,"D[%d,%d]/S[%d,%d] = {%.10f, %.10f, %.10f}\n",i,j,i,j,d_x/S,d_y/S,d_z/S);
    fprintf(target,"_______________________________________________________________\n\n");
//     fflush (target);
    
    if(grab_HS_data){
        int num=i*grabber_n_states+j;
         H_grabbed[num]=H_1el+H_2el+V_nuc*S;
         S_grabbed[num]=S;
        Dx_grabbed[num]=d_x;
        Dy_grabbed[num]=d_y;
        Dz_grabbed[num]=d_z;
    }

    return 0;
}




