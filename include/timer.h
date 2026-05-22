# ifndef TIMER_H
# define TIMER_H
# include <stdio.h>
# include <stdlib.h>
# include <chrono>
# include <vector>

extern std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
extern std::chrono::time_point<std::chrono::high_resolution_clock> t_prev;
extern std::chrono::time_point<std::chrono::high_resolution_clock> t_now;



class mp_job_timer;

class mp_job_timer{
    
public:
    void init(int np);
    void start(int nt);
    void pause(int nt);
    void print(char * name);
    ~mp_job_timer();
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> * t_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> * t_prev;
    std::chrono::time_point<std::chrono::high_resolution_clock> * t_now;
    double * t;
    int n_proc;
    
};


void get_start_time();

void print_current_time(const char * comment);

void enable_print_timers();

void disable_print_timers();

void printf_timer(const char * coment);

void take_timer_data(double * ft,double * ct);

#endif
