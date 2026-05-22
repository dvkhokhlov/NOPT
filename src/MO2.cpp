# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include "keywords.h"
# include "common_vars.h"
# include "inp_out.h"
# include "defaults.h"

int VEC_read(char *inp_name, int start_line, double* vec, int n_ao, int n_mo, char ab){
    
    recursive_file in;
    char * line = new char[BUF_LINE_LENGTH];
    
    in.r_open (inp_name);
    if(in.file[0]==NULL){
        fprintf(out_stream,"Couldn't open '%s'\n",inp_name);
        exit(0);
    }
    
    for(int i=0;i<start_line;i++)in.r_gets(line,BUF_LINE_LENGTH);
    in.r_gets(line,BUF_LINE_LENGTH);
    
    if(ab=='a'){
        while((key_word_comp(line,MO_text_group_start)==0)&&(!in.r_eof())) in.r_gets(line,BUF_LINE_LENGTH);
    }
    else if(ab=='b'){
        while((key_word_comp(line,MO_group_b_start)==0)&&(!in.r_eof())) in.r_gets(line,BUF_LINE_LENGTH);
    }
    else{
        printf("wrong identificator a or b in VEC_read\n");
        exit(1);
    }
    
    if(strstr(line,"[MO]")){
        delete[] line;
        in.r_close();
        return -18;
    }
    
    
    int length=15;
    int first_pos=5;
    
    if(key_word_comp(line,wide)){
        length=25;
        first_pos=6;
        printf("reading VEC-group in wide format\n\n");
//         getchar();
    }
    
    int tmp_mo_num=0;
    int coef_num_row=0;
    
    in.r_gets(line,BUF_LINE_LENGTH);
    while((key_word_comp(line,MO_group_end)==0)&&(!in.r_eof()))
    {
      int tmpint=5;
//       printf("%s\n",line);
      
      if(n_ao<5) tmpint=n_ao;
      
      if(n_ao-coef_num_row*5<5) 
          {
            tmpint=n_ao-coef_num_row*5;
          }
      
      if(n_ao-coef_num_row*5==0) 
      {
          tmp_mo_num++;
          coef_num_row=0;
      }
//       printf("lenth = %d tmpint = %d\n", length, tmpint);
//       getchar();
      char *line2 = new char[length+1];
      
      for(int i=0;i<tmpint;i++)
      {
          for(int j=0;j<length;j++) line2[j]=line[j+length*i+first_pos];
          line2[length] = '\0';
          if(length>15) if(line2[20]=='D') line2[20]='E';
          sscanf(line2,"%lf",&(vec[tmp_mo_num*n_ao+coef_num_row*5+i]));
          
//           if((tmp_mo_num==n_ao-1)&&(coef_num_row*5+i==n_ao-1)) printf("final coef in VEC found\n");
//           printf("%d %d %e\n", tmp_mo_num, coef_num_row*5+i, vec[tmp_mo_num*n_ao+coef_num_row*5+i]); //getchar();
      }
      
      coef_num_row++;
//       printf("%d,%d\n%s",coef_num_row,tmp_mo_num,line);
//       getchar();
      if(n_ao-coef_num_row*5<=0) 
          {
            tmp_mo_num++;
            coef_num_row=0;
          }
      
      in.r_gets(line,BUF_LINE_LENGTH);
      delete[] line2;
    }
//     printf("!!!!!!!!!!!%d\n",n_ao);
//     getchar();
    
//     if(tmp_mo_num!=n_ao) printf("WARNING not all orbitals are given in $VEC group (%d of %d)\n",tmp_mo_num,n_ao);
    
    delete[] line;
    
    in.r_close();
    
    return tmp_mo_num;
}

int energy_read(char *inp_name, int start_line, double* en, int n_mo, char ab){
    
  recursive_file in;
  char line[BUF_LINE_LENGTH];
  
  in.r_open (inp_name);
  for(int i=0;i<start_line;i++)in.r_gets(line,BUF_LINE_LENGTH);
  in.r_gets(line,BUF_LINE_LENGTH);
  if(ab=='a'){
      while((key_word_comp(line,energy_group_start)==0)&&(!in.r_eof())) in.r_gets(line,BUF_LINE_LENGTH);
  }
  else if(ab=='b'){
      while((key_word_comp(line,energy_group_b_start)==0)&&(!in.r_eof())) in.r_gets(line,BUF_LINE_LENGTH);
  }
  else{
      printf("wrong identificator a or b in energy_read\n");
      exit(1);
  }
  
  if(in.r_eof()){
//       printf("WARNING: no orbital energy of type %c in file %s\n\n",ab,inp_name);
      for(int i=0; i<n_mo; i++) en[i]=0.0;
      in.r_close();
      return 0;
  }
  
  int tmp_mo_num=0;
  in.r_gets(line,BUF_LINE_LENGTH);
  while((key_word_comp(line,energy_group_end)==0)&&(!in.r_eof()))
  {
    int tmpint=0;
    sscanf(line,"%lf",&(en[tmp_mo_num]));

    tmp_mo_num++;
    in.r_gets(line,BUF_LINE_LENGTH);
    
  }
//   printf("!!!!!!!!!!!%d\n",n_ao);
//   getchar();
  if(tmp_mo_num!=n_mo) printf("WARNING not all orbitals are given in $energy1 group (%d of %d)\n",tmp_mo_num,n_mo);
  in.r_close();

  return 0;
}
