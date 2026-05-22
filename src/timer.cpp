# include <stdio.h>
# include <chrono> 
# include "timer.h"
# include "common_vars.h"
int print_timers;
std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
std::chrono::time_point<std::chrono::high_resolution_clock> t_prev;
std::chrono::time_point<std::chrono::high_resolution_clock> t_now;

void enable_print_timers(){
    print_timers=1;
}

void disable_print_timers(){
    print_timers=0;
}


void get_start_time(){
    
    t_start = std::chrono::high_resolution_clock::now();
    t_prev = t_start;
    enable_print_timers();
    
}

void print_current_time(const char * comment){
    
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    
    char buffer[100];
    struct tm* tm_info;

    tm_info = localtime(&now_c);

    strftime(buffer, 100, "%H:%M (%Z), %d %B %Y (%A)", tm_info);
    
    
    fprintf(out_stream,"%s%s.\n",comment, buffer);
    
}


void printf_timer(const char * coment){
    
    if(print_timers){
        t_now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_elapsed = t_now - t_prev;
        fprintf(out_stream,"%s:\n", coment);
        fprintf(out_stream,"Time to step - %10.5lf\n", time_elapsed.count());
        time_elapsed = t_now - t_start;
        fprintf(out_stream,"Full time    - %10.5lf\n", time_elapsed.count());
        t_prev=t_now;
        fflush(stdout);
    }
}

void take_timer_data(double * ft,double * ct){

    t_now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = t_now - t_prev;
    
    *ct = time_elapsed.count();
    time_elapsed = t_now - t_start;
    *ft = time_elapsed.count();
    t_prev=t_now;
    
    
}

void mp_job_timer::init(int np){
    
    t_start = new std::chrono::time_point<std::chrono::high_resolution_clock>[np];
    t_prev  = new std::chrono::time_point<std::chrono::high_resolution_clock>[np];
    t_now   = new std::chrono::time_point<std::chrono::high_resolution_clock>[np];
    t = new double[np];
    for(int i=0;i<np;i++) t[i]=0;
    n_proc = np;
}

void mp_job_timer::start(int nt){
    if(nt>=n_proc)fprintf(out_stream,"ERROR IN mp_job_timer.start: process number (%d) is more then max_n_proc set to timer(%d)\n",nt,n_proc);
    t_prev[nt] = std::chrono::high_resolution_clock::now();
}

void mp_job_timer::pause(int nt){
    if(nt>=n_proc)fprintf(out_stream,"ERROR IN mp_job_timer.pause: process number (%d) is more then max_n_proc set to timer(%d)\n",nt,n_proc);
    t_now[nt] = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_elapsed = t_now[nt] - t_prev[nt];
    t[nt]+=time_elapsed.count();
    
}

void mp_job_timer::print(char * name){
    fprintf(out_stream,"Time spent for \"%s\":\n",name);
    fprintf(out_stream,"thread:");
    for(int i=0;i<n_proc;i++){
        if(i<10)fprintf(out_stream,"        %d         |",i);
        else    fprintf(out_stream,"        %d        |",i);
        fprintf(out_stream,"\n");
    }
    fprintf(out_stream,"time  :");
    for(int i=0;i<n_proc;i++){
        fprintf(out_stream," %10.5lf |");
        fprintf(out_stream,"\n");
    }
}
    
mp_job_timer::~mp_job_timer(){
    
    delete[] t_start;
    delete[] t_prev;
    delete[] t_now;
    delete[] t;
}

