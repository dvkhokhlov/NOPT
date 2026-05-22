#include "inp_par_read.h"

int D5=1;
int RI=0;
int RI_flag=-1;
int SO=0;
int nQM=1;
int engines_initialized=0;;
int num_threads=1;
int IS_SYM;
int LINEAR=0;
int write_orbs=1;
int write_ci=0;

// #define NOPT_PRINT_STD_ONLY
FILE * out_stream =stdout;
int os_auto_close=0;
int set_out_stream_by_name(const char * name,const char *mode){

#ifndef NOPT_PRINT_STD_ONLY    
    if(mode[0]!='a')if(mode[0]!='w'){
        fprintf(stdout,"ERROR: unknown parameter \"%c\" (common_vars.cpp)\n",mode[0]);
        fprintf(stdout,"       set_out_stream_by_name can use only \"a\" or \"w\"\n" );
        exit(0);
    }
    
    if(os_auto_close)fclose(out_stream);
    printf("NOPT out_stream: redirected to %s\n", name);
    out_stream = fopen(name,mode);
    
    os_auto_close=1;
#endif
    
    return 0;
}

int set_out_stream_by_FILE(FILE * file, const char * name){
    
#ifndef NOPT_PRINT_STD_ONLY    
    if(os_auto_close)fclose(out_stream);
    
    out_stream = file;
    printf("NOPT out_stream: redirected to previously opened file - %s\n", name);
    os_auto_close=0;
#endif
    
    return 0;
}

int set_out_stream_stdout(){
    
#ifndef NOPT_PRINT_STD_ONLY    
    if(os_auto_close)fclose(out_stream);
    
    out_stream = stdout;
    printf("NOPT out_stream: redirected to stdout\n");
    os_auto_close=0;
#endif
    
    return 0;
}

int set_out_stream_stderr(){
    
#ifndef NOPT_PRINT_STD_ONLY    
    if(os_auto_close)fclose(out_stream);
    
    out_stream = stderr;
    printf("NOPT out_stream: redirected to stderr\n");
    os_auto_close=0;
#endif
    
    return 0;
}

double * RI_AO_M=NULL;

backup_par backup;

int * nopt_parallel_index;

int nopt_parallel_init(int n){
    
    num_threads = n;
    nopt_parallel_index=new int[num_threads+1];
    
    return 0;
}

int nopt_parallel_finalize(){
    
    delete[] nopt_parallel_index;
    
    return 0;
}

