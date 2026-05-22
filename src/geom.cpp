//user
# include "molecule.h"
# include "common_vars.h"
# include "geom.h"
# include "etc.h"
# include "defaults.h"


int shift_molecule(molecule * A, double * dir, double l){
    
    double sc =l/sqrt(dir[0]*dir[0]+
                      dir[1]*dir[1]+
                      dir[2]*dir[2]);
    
    double x = sc*dir[0];
    double y = sc*dir[1];
    double z = sc*dir[2];
    
    for(int i=0;i<A->n_atoms;i++){
        A->atom_coord[3*i+0]+=x*A->au_coef;
        A->atom_coord[3*i+1]+=y*A->au_coef;
        A->atom_coord[3*i+2]+=z*A->au_coef;
    }
    
    for(auto &s: A->s){
        s.O[0]+=x*A->au_coef;
        s.O[1]+=y*A->au_coef;
        s.O[2]+=z*A->au_coef;
    }
    for(auto &s: A->ri_s){
        s.O[0]+=x*A->au_coef;
        s.O[1]+=y*A->au_coef;
        s.O[2]+=z*A->au_coef;
    }
    
    
    return 0;
}


coord_list C;

coord_list::coord_list(){
    
    value.resize(0);
    refs.resize(0);
    name.resize(0);
    
    
}

coord_list::~coord_list(){
    

}


int read_coord(char * l, double * v, int p, int t){
    
    int s=0;
    while(l[s]==' ')s++;
    
    if(isdigit(l[s])||(l[s]=='-')){
        sscanf(l+s,"%lf", v+p);
        while(isdigit(l[s])||(l[s]=='.'))s++;
        if(t)v[p]=v[p]*M_PI/180.0;
        
        return s;
    }
    
    if(isalpha(l[s])){
        std::vector<char> c;
        c.resize(0);
        while (isalpha(l[s])||isdigit(l[s])){
            c.push_back(l[s]);
            s++;
        }
        c.push_back('\0');
        
        for(int i=0;i<C.name.size();i++){
            if(c==C.name[i]){
                if(t==C.type[i]){
                    C.refs[i].push_back(p);
                    return s;
                }
                else{
                    printf("ERROR: variable %s is used for both distance and angle\n",c.data());
                    exit(0);
                }
            }
        }
        C.name.push_back(c);
        C.type.push_back(t);
        std::vector<int> pos;
        pos.resize(1);
        pos[0]=p;
        C.refs.push_back(pos);
        
        return s;
    }
    
    printf("ERROR: read_coord couldn't read number\n");
    
    exit(0);
    
}


int read_cl(recursive_file * in){
    
    char * line = new char[BUF_LINE_LENGTH];
    
    char * l;
    C.value.resize(C.name.size());
    in->r_gets(line,BUF_LINE_LENGTH);
    while(in->r_eof()!=0){
        for(int i=0;i<C.name.size();i++){

            l=strstr(line,C.name[i].data());
            
            if(l!=nullptr) {
                l+=C.name[i].size();
                sscanf(l, "%lf",&(C.value[i]));
                if(C.type[i])
                    C.value[i]=C.value[i]*M_PI/180.0;
            }
        }
        in->r_gets(line,BUF_LINE_LENGTH);
    }
    
    
    
    delete[] line;
        
    return 0;
}

int dist_scale(double c){
    
    for(int i=0;i<C.value.size();i++)
        if(C.type[i]==0)
            C.value[i]=C.value[i]*c;
    
    return 0;
}


int coord_sync(double * v){
    
    for(int i=0;i<C.value.size();i++){
        for(auto &p: C.refs[i])
            v[p]=C.value[i];
    }
    
    
    return 0;
}

int print_coord_list(){
    
    double v;
    for(int i=0;i<C.value.size();i++){
        if(C.type[i])v=C.value[i]/M_PI*180;
        else         v=C.value[i]/1.889725989;
        fprintf(out_stream,"%s = %f\n",C.name[i].data(),v);
    }
    
    getchar();
    
    return 0;
    
}

int fprint_coord_list(const char * name){
    
    FILE * out_file = fopen(name,"w");
    FILE * backup;
    backup =out_stream;
    out_stream = out_file;
    print_coord_list();
    out_stream = backup;
    fclose(out_file);
    
    return 0;
}


int set_std_delta(){
    
    
    int n = C.value.size();
    C.value_delta.resize(n);
    
    double dR=0.002;
    double dA=0.002;//*M_PI/180.0;
    
    for(int i=0;i<n;i++){
        if(C.type[i])C.value_delta[i]=dA;
        else         C.value_delta[i]=dR;
        
    }
    
    return 0;
    
}
